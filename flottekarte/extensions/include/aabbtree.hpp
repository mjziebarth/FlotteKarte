/*
 * Axis-aligned bounding box tree.
 *
 * Authors: Malte J. Ziebarth (malte.ziebarth@tum.de)
 *
 * Copyright (C) 2024 Technische Universität München
 *
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
 * the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 */

#ifndef FLOTTEKARTE_AABBTREE_HPP
#define FLOTTEKARTE_AABBTREE_HPP

#include <limits>
#include <iostream>
#include <thread>
#include <chrono>
#include <../include/types.hpp>
#include <../include/geometry.hpp>

/*
 * This contains some more detail that would clog this header:
 */
#include <../include/backbone/aabbbranch.hpp>

namespace flottekarte {


template<typename index_t=size_t>
class AABBTree
{
public:

    AABBTree()
    {};


    /*
     * Adding segments:
     * ----------------
     */
    void add_segment(const segment_t& s);

    void add_segment(const xy_t& p0, const xy_t& p1)
    {
        add_segment(segment_t(p0, p1));
    }

    /*
     * Query whether a geometry intersects.
     * ------------------------------------
     *
     * Intersection with segments:
     */
    bool intersects(const segment_t& s) const;

    bool intersects(const xy_t& p0, const xy_t& p1) const
    {
        return intersects(segment_t(p0, p1));
    }

    /*
     * Intersection with circles:
     */
    bool intersects(const circle_t& c) const;

    /*
     * Iterate all elements of an intersection:
     */
    aabb::Intersection<segment_t, index_t>
    intersection(const segment_t& s) const;

    aabb::Intersection<circle_t, index_t>
    intersection(const circle_t& c) const;


    /*
     * Debug functionality
     * -------------------
     */
    bool internally_valid() const;
    void debug_print() const;

private:
    typedef aabb::branch_t<index_t> branch_t;
    typedef aabb::leaf_t<index_t> leaf_t;

    constexpr static index_t root_jid = aabb::BRANCH0<index_t>;

    /* We use 'size_t' to index_t both in branch_t and leaf_t vectors.
     * This value is used to separate the indices in both vectors: */

    std::vector<branch_t> branches;
    std::vector<leaf_t> leaves;

    index_t find_insert_sibling(const bbox_t& bbox) const;

    bbox_t get_bbox(index_t jid) const;

    const branch_t& root() const;

    const branch_t& branch(index_t jid) const;

    branch_t& branch(index_t jid);

    void debug_print_impl() const;
    bool internally_valid_impl() const;



    /*
     * Ensure that the intersections are ranges:
     */
    static_assert(std::ranges::input_range<
                     aabb::Intersection<segment_t,index_t>
                  >,
                  "Segment intersection not a range");
    static_assert(std::ranges::input_range<
                     aabb::Intersection<segment_t,index_t>
                  >,
                  "Circle intersection not a range");
};



/*****************************************************************************
 *
 *                             Implementation
 *
 *****************************************************************************/


template<typename index_t>
bool AABBTree<index_t>::intersects(const segment_t& s) const
{

    /* Early exits: */
    if (leaves.size() == 0){
        return false;
    } else if (leaves.size() == 1){
        return flottekarte::intersects(leaves[0].segment, s);
    }


    /* Otherwise traverse the tree: */
    typedef aabb::Intersection<segment_t, index_t> Intersection;
    return !Intersection(branches, leaves, s).empty();
}


template<typename index_t>
bool AABBTree<index_t>::intersects(const circle_t& c) const
{
    /* Early exits: */
    if (leaves.size() == 0){
        return false;
    } else if (leaves.size() == 1){
        return flottekarte::intersects(leaves[0].segment, c);
    }

    /* Otherwise traverse the tree:: */
    typedef aabb::Intersection<circle_t, index_t> Intersection;
    return !Intersection(branches, leaves, c).empty();
}


template<typename index_t>
aabb::Intersection<segment_t, index_t>
AABBTree<index_t>::intersection(const segment_t& s) const
{
    return aabb::Intersection<segment_t, index_t>(branches, leaves, s);
}


template<typename index_t>
aabb::Intersection<circle_t, index_t>
AABBTree<index_t>::intersection(const circle_t& c) const
{
    return aabb::Intersection<circle_t, index_t>(branches, leaves, c);
}


template<typename index_t>
void AABBTree<index_t>::add_segment(const segment_t& s)
{
    /* Ensure that the current index_t can represent all segments: */
    index_t lid = leaves.size();
    if (lid >= aabb::NONE<index_t>)
        throw std::runtime_error("Maximum number of segments in AABBTree "
            "exceeded."
        );

    /* Add the leaf to the vector: */
    leaves.emplace_back(s);
    bbox_t bbox(leaves.back().bbox());


    /* Some early exits: */
    if (lid == 0){
        /* First leaf. */
        return;
    } else if (lid == 1){
        /* Second leaf. Create the first branch: */
        branches.emplace_back(
            bbox_t(leaves[0].bbox(), bbox),
            aabb::NONE<index_t>, 0, 1
        );
        leaves[0].parent = aabb::bid2jid<index_t>(0);
        leaves[1].parent = aabb::bid2jid<index_t>(0);
        return;
    }



    /* Use heuristic to find good sibling: */
    index_t sibling_jid = find_insert_sibling(bbox);
    if (!aabb::is_leaf(sibling_jid))
        throw std::runtime_error("find_insert_sibling did not return leaf.");
    if (sibling_jid >= leaves.size())
        throw std::runtime_error("AABB find_insert_sibling returns corrupt "
            "sibling jid (too large)."
        );

    /* New branch ID, and check whether sufficient indices are available:
     * Note: Not sure if this is even possible or whether the additional
     * leaf will always trigger index_t exceedence first.
     */
    index_t new_bid = branches.size();
    constexpr index_t max_index_t = std::numeric_limits<index_t>::max();
    if (new_bid >= max_index_t - aabb::BRANCH0<index_t>)
        throw std::runtime_error("Maximum number of branches in AABBTree "
            "exceeded."
        );

    /* Add a new branch at the position of the sibling, including the
     * newly added leaf and the sibling. */
    leaf_t& sibling = leaves.at(sibling_jid);
    index_t grandparent_jid = sibling.parent;

    bbox_t parent_bbox(bbox, sibling.bbox());
    branches.emplace_back(parent_bbox, grandparent_jid, sibling_jid, lid);
    index_t new_jid = aabb::bid2jid(new_bid);

    /* The newly added branch is parent to both its leaves: */
    sibling.parent = new_jid;
    leaves.back().parent = new_jid;

    /* Need to adjust grandparent's child IDs and (potentially) bbox: */
    branch_t& grandparent = branch(grandparent_jid);
    if (grandparent.c0 == sibling_jid)
        grandparent.c0 = new_jid;
    else if (grandparent.c1 == sibling_jid)
        grandparent.c1 = new_jid;
    else
        throw std::runtime_error("AABBTree data corrupted: child not in parent "
            "ids!"
        );


    grandparent.bbox = bbox_t(parent_bbox, grandparent.bbox);

    /* Bottom-Up increase of the depth and adjustment of the bbox: */
    branch_t* ancestor = &grandparent;
    bbox_t last_bbox = grandparent.bbox;
    while (ancestor->parent != aabb::NONE<index_t>)
    {
        ancestor = &branch(ancestor->parent);
        if (contains(ancestor->bbox, last_bbox))
            break;
        last_bbox = bbox_t(ancestor->bbox, last_bbox);
        ancestor->bbox = last_bbox;
    }
}


template<typename index_t>
index_t AABBTree<index_t>::find_insert_sibling(const bbox_t& bbox) const
{
    if (branches.empty())
        throw std::runtime_error("Trying to find an insert sibling when there "
            "is no branch yet."
        );
    /* This boolean tries to alternate ties: */
    bool choose_left = (leaves.size() % 2) == 0;

    /* Start from the root branch and descend until we have found a leaf.
     * To find a suitable node, we use the heuristic described in
     * https://www.azurefromthetrenches.com/introductory-guide-to-aabb-tree-collision-detection/
     * that is, choosing the node which minimizes the joint area with 'bbox'.
     */
    index_t jid = root_jid;
    while (aabb::is_branch(jid))
    {
        /* Get the current branch: */
        const branch_t& b = branch(jid);

        /* Check the two bounding boxes of the children: */
        bbox_t b0(get_bbox(b.c0));
        bbox_t b1(get_bbox(b.c1));

        /* Compute the joint areas: */
        double A0 = bbox_t(b0, bbox).area();
        double A1 = bbox_t(b1, bbox).area();

        if (A0 < A1)
            jid = b.c0;
        else if (A1 < A0)
            jid = b.c1;
        else if (choose_left) {
            choose_left = false;
            jid = b.c0;
        } else {
            choose_left = true;
            jid = b.c1;
        }
    }

    /* We find ourselves on a leaf now. */
    return jid;
}


template<typename index_t>
bbox_t AABBTree<index_t>::get_bbox(index_t jid) const
{
    if (aabb::is_branch(jid))
        return branch(jid).bbox;
    else if (aabb::is_leaf(jid))
        return leaves.at(jid).bbox();

    throw std::runtime_error("Trying to get bbox of invalid joint ID.");

}


template<typename index_t>
const aabb::branch_t<index_t>& AABBTree<index_t>::root() const
{
    if (branches.empty())
        throw std::runtime_error("Trying to access root branch when there are "
            "no branches yet."
        );

    return branches.front();
}


template<typename index_t>
const aabb::branch_t<index_t>& AABBTree<index_t>::branch(index_t jid) const
{
    return branches.at(jid - aabb::BRANCH0<index_t>);
}


template<typename index_t>
aabb::branch_t<index_t>& AABBTree<index_t>::branch(index_t jid)
{
    return branches.at(jid - aabb::BRANCH0<index_t>);
}


} // namespace flottekarte



/*
 * Implementation:
 */
namespace flottekarte {

template<>
bool AABBTree<size_t>::internally_valid() const;

template<>
void AABBTree<size_t>::debug_print() const;

}



#endif
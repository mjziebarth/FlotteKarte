/*
 * Axis-aligned bounding box tree: backbone detail
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

#ifndef FLOTTEKARTE_BACKBONE_AABBTREE_HPP
#define FLOTTEKARTE_BACKBONE_AABBTREE_HPP

#include <../include/types.hpp>
#include <../include/geometry.hpp>

#include <ranges>

namespace flottekarte {

namespace aabb {

template<typename index_t>
constexpr static index_t NONE = std::numeric_limits<index_t>::max()
                                / static_cast<index_t>(2);

template<typename index_t>
constexpr static index_t BRANCH0 = NONE<index_t> + static_cast<index_t>(1);

template<typename index_t>
struct branch_t;

/*
 * leaf_t
 * ------
 */

template<typename index_t>
struct leaf_t {
    typedef index_t index_type;

    segment_t segment;
    index_t parent;

    leaf_t() = default;
    leaf_t(const segment_t& s) : segment(s), parent(NONE<index_t>)
    {}

    bbox_t bbox() const
    {
        return bbox_t(segment);
    }
};

template<typename index_t>
bool contains(const bbox_t& b, const leaf_t<index_t>& l)
{
    return contains(b, l.bbox());
}



/*
 * branch_t
 * --------
 */

template<typename index_t>
struct branch_t {
    typedef index_t index_type;

    bbox_t bbox;
    index_t parent;
    index_t c0;
    index_t c1;

    branch_t(const bbox_t& bbox, index_t parent, index_t c0,
             index_t c1)
       : bbox(bbox), parent(parent), c0(c0), c1(c1)
    {}


    static branch_t root()
    {
        return branch_t();
    }

private:
    branch_t() : bbox(xy_t(0,0), xy_t(0,0)), parent(NONE<index_t>)
    {}
};


/*
 * handle jids (=joint ID):
 */
template<typename index_t>
constexpr static bool is_branch(index_t jid){
    return jid >= BRANCH0<index_t>;
}

template<typename index_t>
constexpr static bool is_leaf(index_t jid){
    return jid < NONE<index_t>;
}


template<typename index_t>
constexpr static index_t bid2jid(index_t bid){
    return bid + BRANCH0<index_t>;
}


template<typename index_t>
constexpr static index_t jid2bid(index_t jid){
    return jid - BRANCH0<index_t>;
}




/*
 * Sentry for the IntersectionIterator:
 */
template<typename T, typename index_t>
class IntersectionIterator;

struct IntersectionIteratorEnd
{
    template<typename T, typename index_t>
    bool operator==(const IntersectionIterator<T,index_t>& it) const
    {
        return it.is_done;
    }
};

/*
 * IntersectionIterator
 * --------------------
 */

template<typename T, typename index_t>
class IntersectionIterator
{
    friend IntersectionIteratorEnd;
public:
    /* Iterator traits: */
    typedef segment_t value_type;
    typedef const segment_t* pointer;
    typedef const segment_t& reference;
    typedef std::ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;

    typedef flottekarte::aabb::branch_t<index_t> branch_t;
    typedef flottekarte::aabb::leaf_t<index_t> leaf_t;


    IntersectionIterator(
        const std::vector<branch_t>& branches,
        const std::vector<leaf_t>& leaves,
        const T& _geometry
    )
       : branches(branches), leaves(leaves), is_done(false),
         current(BRANCH0<index_t>), geometry(_geometry),
         geometry_bbox(_geometry)
    {
        /* Special cases: */
        if (leaves.empty()){
            /* No leaf - no intersection. */
            is_done = true;
            return;
        }
        if (!intersects(branches[0].bbox, geometry_bbox)){
            /* There is no intersection! */
            std::cout << "globally no intersection.\n";
            is_done = true;
            return;
        }
        if (leaves.size() == 1){
            /* Just one leaf but no branch structure so far.
             * Indicate that we cannot advance branches: */
            branch_stack.push_back(2);

            /* Dereference to the single leaf: */
            current = 0;

            return;
        }

        /* "Standard" case: */

        /* Start with the first branch: */
        branch_stack.push_back(0);

        /* Advance till we found the first intersection: */
        operator++();
    }

    IntersectionIterator(const IntersectionIterator& other) = default;
    IntersectionIterator(IntersectionIterator&& other) = default;

    ~IntersectionIterator() noexcept(true) {};


    static IntersectionIterator end(
        const std::vector<branch_t>& branches,
        const std::vector<leaf_t>& leaves
    )
    {
        return IntersectionIterator(branches, leaves);
    }

    IntersectionIterator& operator=(const IntersectionIterator& other) = default;
    IntersectionIterator& operator=(IntersectionIterator&& other) = default;

    bool operator==(const IntersectionIterator& other) const
    {
        /* Ensure that both iterators iterate the same tree: */
        if (&branches.get() != &other.branches.get())
            return false;

        /* The end iterator is special and is tagged by 'is_done'. */
        if (is_done && other.is_done)
            return true;

        /* Else compare: */
        return (geometry == other.geometry) && (current == other.current);
    }


    bool operator==(const IntersectionIteratorEnd&) const
    {
        return is_done;
    }


    IntersectionIterator& operator++()
    {
        if (current == NONE<index_t> || is_done)
            throw std::runtime_error("Trying to advance 'NONE' "
                "IntersectionIterator::operator++()"
            );

        /* Iterate until we have exhausted both branches of the root node: */
        while (current != BRANCH0<index_t> || branch_stack.front() != 2)
        {
            /* Here we perform one advancement step. */
            if (is_leaf(current)){
                /* We need to go up. */
                branch_stack.pop_back();
                current = leaves.get()[current].parent;
                continue;
            } else {
                index_t next_jid;
                const branch_t& b(branches.get()[jid2bid(current)]);
                if (branch_stack.back() == 0){
                    /* We want to descend to the left branch. */
                    next_jid = b.c0;

                } else if (branch_stack.back() == 1){
                    /* We want to descend to the right branch. */
                    next_jid = b.c1;

                } else {
                    /* Ascend. */
                    branch_stack.pop_back();
                    current = b.parent;
                    continue;
                }
                /* Advance in this layer */
                ++branch_stack.back();

                /* Check whether the proposed branch intersects: */
                if (is_branch(next_jid)){
                    const branch_t& bnext(
                        branches.get()[jid2bid(next_jid)]
                    );
                    if (!intersects(bnext.bbox, geometry_bbox))
                        /* No intersection - already clear from the bbox.
                         * Continue! */
                        continue;

                } else {
                    const leaf_t& lnext(leaves.get()[next_jid]);
                    if (!intersects(lnext.bbox(), geometry_bbox))
                        /* No intersection - already clear from the bbox.
                         * Continue! */
                        continue;
                }

                /* Add the new layer: */
                branch_stack.push_back(0);
                current = next_jid;
            }

            /* Now check if the new node is a leaf.
             * If so, check whether it intersects: */
            if (is_leaf(current)){
                const leaf_t& l(leaves.get()[current]);
                if (intersects(l.segment, geometry))
                {
                    /* Found an intersection! */
                    return *this;
                }
            }
        }
        /* If we arrived here, we are done! */
        is_done = true;
        current = NONE<index_t>;

        return *this;
    }

    IntersectionIterator& operator++(int)
    {
        IntersectionIterator copy = *this;
        operator++();
        return copy;
    }

    const segment_t& operator*() const
    {
        if (is_done)
            throw std::runtime_error("Trying to dereference "
                "IntersectionIterator at end."
            );
        if (!aabb::is_leaf(current))
            throw std::runtime_error("Internal state in IntersectionIterator "
                "corrupted: dereferencing at non-leaf index."
            );
        return leaves.at(current).leaf;
    }

    index_t leaf_id() const
    {
        return current;
    }

private:
    std::reference_wrapper<const std::vector<branch_t>> branches;
    std::reference_wrapper<const std::vector<leaf_t>> leaves;
    bool is_done;
    index_t current;
    std::vector<uint8_t> branch_stack;
    T geometry;
    bbox_t geometry_bbox;

    IntersectionIterator(
        const std::vector<branch_t>& branches,
        const std::vector<leaf_t>& leaves
    )
       : branches(branches), leaves(leaves), is_done(true),
         current(NONE<index_t>), geometry_bbox(xy_t(0,0), xy_t(0,0))
    {}

};



/*
 * The generator for the previous iterator:
 */
template<typename T, typename index_t>
class Intersection
   : public std::ranges::view_interface<IntersectionIterator<T,index_t>>
{
public:
    typedef IntersectionIterator<T, index_t> iterator;

    typedef flottekarte::aabb::branch_t<index_t> branch_t;
    typedef flottekarte::aabb::leaf_t<index_t> leaf_t;


    Intersection(
        const std::vector<branch_t>& branches,
        const std::vector<leaf_t>& leaves,
        const T& geometry
    )
       : _begin(branches, leaves, geometry)
    {}

    IntersectionIterator<T,index_t> begin() const
    {
        return _begin;
    }

    IntersectionIteratorEnd end() const
    {
        return IntersectionIteratorEnd();
    }

    bool empty() const
    {
        /* This is not thread safe... except that it probably is?
         * After all, _empty will just be overridden by the same boolean
         * variable and all that may happen is redundant computation.
         * That is, if there is no race condition problem in _empty.emplace
         * other than overwriting previous values. */
        if (!_empty.has_value())
            _empty.emplace(begin() == end());
        return *_empty;
    }

private:
    /* Cache the begin iterator to achieve amortized constant access: */
    IntersectionIterator<T, index_t> _begin;

    mutable std::optional<bool> _empty;

    /*
     * Ranges concepts.
     */
    static_assert(std::destructible<iterator>, "Not destructible.");
    static_assert(std::constructible_from<iterator, iterator>,
                  "Not constructible from.");
    static_assert(std::move_constructible<iterator>, "Not move-constructible.");
    static_assert(std::assignable_from<iterator&, iterator>,
                  "Not assignable from.");
    static_assert(std::swappable<iterator>, "Not swappable.");
    static_assert(std::movable<iterator>, "Not movable.");
    static_assert(std::weakly_incrementable<iterator>,
                  "Not weakly incrementable.");

    static_assert(std::input_or_output_iterator<iterator>,
                  "Not an input_or_output_iterator.");
    static_assert(std::input_iterator<iterator>,
                  "Not an input_iterator.");

    static_assert(std::sentinel_for<IntersectionIteratorEnd, iterator>,
                  "end() not a sentinel for iterator.");
};

} // namespace aabb

} // namespace flottekarte

#endif
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

#define FLOTTEKARTE_AABBTREE_CPP

#include <../include/aabbtree.hpp>

#include <iostream>

namespace flottekarte {





/*
 * AABBTree
 * --------
 */



/*
 * A debugging test.
 */
template<typename index_t>
bool AABBTree<index_t>::internally_valid_impl() const
{
    if (branches.empty()){
        if (leaves.empty())
            return true;
        return leaves.front().parent == aabb::NONE<index_t>;
    }
    auto check_child_links = [this](index_t jid) -> bool
    {
        if (aabb::is_leaf(jid))
            return true;
        const branch_t& b(branches[aabb::jid2bid(jid)]);
        for (index_t jid_c : {b.c0, b.c1}){
            if (aabb::is_leaf(jid_c)){
                if (leaves[jid_c].parent != jid){
                    std::cerr << "leaf.parent != jid\n";
                    std::cerr << "leaf:        " << jid_c << "\n";
                    std::cerr << "leaf.parent: "
                        << aabb::jid2bid(leaves[jid_c].parent) << "\n";
                    std::cerr << "bid:         "
                        << aabb::jid2bid(jid) << "\n";
                    return false;
                }
            } else {
                if (branches[aabb::jid2bid(jid_c)].parent != jid){
                    std::cerr << "branch.parent != jid\n";
                    std::cerr << "branch.parent: "
                        << branches[aabb::jid2bid(jid_c)].parent << "\n";
                    std::cerr << "jid:         " << jid << "\n";
                    return false;
                }
            }
        }
        return true;
    };

    auto check_bbox = [this](index_t jid) -> bool
    {
        if (aabb::is_leaf(jid))
            return true;
        const branch_t& b(branches[aabb::jid2bid(jid)]);
        for (index_t jid_c : {b.c0, b.c1}){
            if (aabb::is_leaf(jid_c)){
                if (!contains(b.bbox, leaves[jid_c].bbox())){
                    std::cerr << "parent bbox does not contain leaf bbox\n";
                    std::cerr << "leaf:        " << jid_c << "\n";
                    std::cerr << "bid:         " << aabb::jid2bid(jid) << "\n";
                    bbox_t lbb(leaves[jid_c].bbox());
                    std::cout << std::setprecision(18);
                    std::cout << "bbox[leaf]:  (" << lbb.xmin << ", "
                        << lbb.xmax << ", " << lbb.ymin << ", " << lbb.ymax
                        << ")\n" << std::flush;
                    std::cout << "bbox[parent]:(" << b.bbox.xmin << ", "
                        << b.bbox.xmax << ", " << b.bbox.ymin << ", "
                        << b.bbox.ymax << ")\n" << std::flush;
                    return false;
                }
            } else {
                if (!contains(b.bbox, branches[aabb::jid2bid(jid_c)].bbox)){
                    std::cerr << "parent bbox does not contain child bbox\n";
                    std::cerr << "bid:         " << aabb::jid2bid(jid) << "\n";
                    return false;
                }
            }
        }
        return true;
    };

    for (index_t bid=1; bid < branches.size(); ++bid)
    {
        if (!check_child_links(aabb::bid2jid(bid))){
            std::cerr << "BRANCH [" << bid << "] IS NOT VALID!\n";
            return false;
        }
        if (!check_bbox(aabb::bid2jid(bid))){
            std::cerr << "BRANCH [" << bid << "] BBOX NOT CONTAINED!\n";
            return false;
        }
    }

    return true;
}


/*
 * Explicit specializations:
 */
template<>
bool AABBTree<size_t>::internally_valid() const
{
    return internally_valid_impl();
}


/*
 * Debugging info to reconstruct the tree.
 */
template<typename index_t>
void AABBTree<index_t>::debug_print_impl() const
{
    std::cout << "AABBTree{\n   branches: [";
    for (index_t i=0; i<branches.size(); ++i){
        const branch_t& b = branches[i];
        std::cout << "(" << i << " : p=" << aabb::jid2bid(b.parent)
            << "; c=[";
        if (aabb::is_leaf(b.c0))
            std::cout << "l" << b.c0 << ",";
        else
            std::cout << "b" << aabb::jid2bid(b.c0) << ",";
        if (aabb::is_leaf(b.c1))
            std::cout << "l" << b.c1 << "]";
        else
            std::cout << "b" << aabb::jid2bid(b.c1) << "]";
        std::cout << "; bbox=[" << b.bbox.xmin << ", " << b.bbox.xmax << ", "
                  << b.bbox.ymin << ", " << b.bbox.ymax << "]), ";
    }
    std::cout << "],\n   leaves: [";

    for (index_t i=0; i<leaves.size(); ++i){
        const leaf_t& l = leaves[i];
        bbox_t bb(l.bbox());
        std::cout << "(" << i << " : p=" << aabb::jid2bid(l.parent)
            << "; bbox=[" << bb.xmin << ", " << bb.xmax << ", "
            << bb.ymin << ", " << bb.ymax << "]),\n";
    }
    std::cout << "]\n}\n" << std::flush;
}


template<>
void AABBTree<size_t>::debug_print() const
{
    debug_print_impl();
}


}
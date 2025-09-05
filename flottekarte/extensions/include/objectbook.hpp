/*
 * Bookkeeping of objects by integer identifiers. This is mostly useful
 * for type-safe, low-overhead object reference communication to the
 * Python world via ctypes.
 *
 * Author: Malte J. Ziebarth (malte.ziebarth@tum.de)
 *
 * Copyright (C) 2025 Technische Universität München
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

#ifndef FLOTTEKARTE_OBJECTBOOK_HPP
#define FLOTTEKARTE_OBJECTBOOK_HPP

#include <variant>
#include <unordered_map>
#include <mutex>
#include <limits>

namespace flottekarte {


/*
 * Transform a set of types to a set of shared_ptr of
 * the same type.
 */
template<typename... Objects>
struct to_shared_ptr_variants;


template<typename O, typename... Objects>
struct to_shared_ptr_variants<O, Objects...>
{
    using type = std::decay_t<
        decltype(
            std::tuple_cat(
                std::declval<std::tuple<std::shared_ptr<O>>>(),
                std::declval<typename to_shared_ptr_variants<Objects...>::type>()
            )
        )
    >;
};

template<typename O>
struct to_shared_ptr_variants<O>
{
    using type = std::tuple<std::shared_ptr<O>>;
};

template<typename... Objects>
using to_shared_ptr_variants_t = typename to_shared_ptr_variants<Objects...>::type;

/*
 * Transform a tuple of types to a variant of the same types.
 */
template<typename... T>
struct tuple_to_variant;

template<typename... T>
struct tuple_to_variant<std::tuple<T...>>
{
    using type = std::variant<T...>;
};

template<typename... T>
using tuple_to_variant_t = typename tuple_to_variant<T...>::type;


/*
 * Ensure that a type is part of a template parameter pack.
 */
template<typename Tq, typename... T>
struct contains_type;

template<typename Tq, typename T0, typename... T>
struct contains_type<Tq, T0, T...>
{
    constexpr static bool value =
        std::is_same_v<T0, Tq>
        || contains_type<Tq, T...>::value;
};

template<typename Tq, typename T0>
struct contains_type<Tq, T0>
{
    constexpr static bool value = std::is_same_v<T0, Tq>;
};

template<typename Tq, typename... T>
constexpr bool contains_type_v = contains_type<Tq, T...>::value;




/*
 * The ObjectBook class
 * ====================
 */
template<typename Int, typename... Objects>
class ObjectBook
{
private:
    constexpr static size_t MAX_INT = std::numeric_limits<Int>::max();
    constexpr static size_t MAX_SIZE_T = std::numeric_limits<size_t>::max();


public:
    using value_type = tuple_to_variant_t<to_shared_ptr_variants_t<Objects...>>;


    std::shared_ptr<const value_type> operator()(Int i) const
    {
    }


    std::shared_ptr<value_type> operator()(Int i)
    {
    }


    template<typename T>
    std::shared_ptr<T> get(Int i)
    {
        auto it = i2o.find(i);
        if (it == i2o.end()){
            std::string msg("No object of ID '");
            msg += std::to_string(i);
            msg += "' exists.";
            throw std::runtime_error(msg);
        }
        return std::visit(
            [i](auto&& ptr) -> std::shared_ptr<T>
            {
                using Tp_in = std::decay_t<decltype(ptr)>;
                using Tp = std::shared_ptr<T>;
                if constexpr (std::is_same_v<Tp_in, Tp>)
                    return ptr;
                std::string msg("Object of ID '");
                msg += std::to_string(i);
                msg += "' is not of the requested type.";
                throw std::runtime_error(msg);
            },
            it->second
        );
    }


    template<typename T>
    std::shared_ptr<const T> get(Int i) const
    {
        auto it = i2o.find(i);
        if (it == i2o.end()){
            std::string msg("No object of ID '");
            msg += std::to_string(i);
            msg += "' exists.";
            throw std::runtime_error(msg);
        }
        return std::visit(
            [i](auto&& ptr) -> std::shared_ptr<const T>
            {
                using Tp_in = std::decay_t<decltype(ptr)>;
                using Tp = std::shared_ptr<T>;
                if constexpr (std::is_same_v<Tp_in, Tp>)
                    return ptr;
                std::string msg("Object of ID '");
                msg += std::to_string(i);
                msg += "' is not of the requested type.";
                throw std::runtime_error(msg);
            },
            it->second
        );
    }


    template<typename Object>
    Int emplace(Object&& o)
    {
        static_assert(contains_type_v<Object, Objects...>);

        /* Lock: */
        std::lock_guard lock(mutex);

        /* First get the ID: */
        size_t oid;
        if (free_ids.empty()){
            oid = i2o.size();
        } else {
            oid = free_ids.back();
        }

        /* Ensure that this index can be expressed in the supplied integer
        * range: */
        if constexpr (MAX_INT < MAX_SIZE_T)
            if (oid > static_cast<size_t>(MAX_INT))
                throw std::runtime_error(
                    "Number of objects exceeds maximum number of ObjectBook "
                    "integer type."
                );

        /* Emplace the object. Here we ensure that the resulting
         * pointer is not nullptr. */
        std::shared_ptr<Object> ptr(std::make_shared<Object>(std::move(o)));
        if (!ptr)
            throw std::runtime_error(
                "Failed to make shared pointer of object."
            );
        i2o.emplace(oid, ptr);

        /* If we got the ID from free_ids, mark that ID as used: */
        if (free_ids.empty())
            free_ids.pop_back();

        return oid;
    }


    void erase(Int oid)
    {
        /* Lock: */
        std::lock_guard lock(mutex);

        /* See if the object ID is valid: */
        auto it = i2o.find(oid);
        if (it == i2o.end()){
            std::string msg("Trying to erase object of ID '");
            msg += std::to_string(oid);
            msg += "' which is not valid.";
            throw std::runtime_error(msg);
        }

        /* Delete the object: */
        i2o.erase(it);

        /* Free the object ID: */
        free_ids.push_back(oid);
    }

private:
    /*
     * Data members:
     */
    std::unordered_map<Int, value_type> i2o;
    std::vector<Int> free_ids;
    std::mutex mutex;


};



}


#endif
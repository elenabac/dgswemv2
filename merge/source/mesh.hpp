#ifndef MESH_HPP
#define MESH_HPP

#include "general_definitions.hpp"
#include "utilities/heterogeneous_containers.hpp"

namespace Geometry {
	template <typename T, typename... Args>
	T create(Args&&... args) {
		return T(std::forward<Args>(args)...);
	}

	//Since elements types already come in a tuple. We can use specialization
	//to get easy access to the parameter packs for the element and edge types.
	template<typename ElementTypeTuple, typename InteriorEdgeTypeTuple, typename BoundaryEdgeTypeTuple>
	class Mesh;

	template<typename... Elements, typename... Interfaces, typename... Boundaries>
	class Mesh<std::tuple<Elements...>, std::tuple<Interfaces...>, std::tuple<Boundaries...>> {
	private:
		using ElementContainer = Utilities::HeterogeneousMap<Elements...>;
		using InterfaceContainer = Utilities::HeterogeneousVector<Interfaces...>;
		using BoundaryContainer = Utilities::HeterogeneousVector<Boundaries...>;

		ElementContainer elements;
		InterfaceContainer interfaces;
		BoundaryContainer boundaries;
	
	public:
		uint GetNumberElements() { return this->elements.size(); }
		uint GetNumberInterfaces() { return this->interfaces.size(); }
		uint GetNumberBoundaries() { return this->boundaries.size(); }
		
		template<typename T, typename... Args>
		void CreateElement(uint ID, Args&&... args) {
			this->elements.template emplace<T>(ID, create<T>(std::forward<Args>(args)...));
		}

		template<typename T, typename... Args>
		void CreateInterface(Args&&... args) {
			this->interfaces.template emplace_back<T>(create<T>(std::forward<Args>(args)...));
		}

		template<typename T, typename... Args>
		void CreateBoundary(Args&&... args) {
			this->boundaries.template emplace_back<T>(create<T>(std::forward<Args>(args)...));
		}

		/*
		template<typename T>
		T* find_element(int id)
		{
			T* ptr = nullptr;
			auto& elements_of_type_T = std::get< util::index<T, TupleType>::value >(_elements.data);
			for (uint i = 0; i < elements_of_type_T.size(); ++i) {
				if (elements_of_type_T[i].get_id() == id) {
					ptr = &(elements_of_type_T[i]);
				}
			}
			return ptr;
		}

		template<typename UnaryFunction>
		void call_for_each_element(UnaryFunction f)
		{
			util::for_each_in_tuple(_elements.data,
				[f](auto& element_vector) {
				std::for_each(element_vector.begin(), element_vector.end(),
					f);
			});
		}

		template<typename UnaryFunction>
		void call_for_each_interior_edge(UnaryFunction f)
		{
			util::for_each_in_tuple(_interior_edges.data,
				[f](auto& edge_vector) {
				std::for_each(edge_vector.begin(), edge_vector.end(),
					f);
			});
		}

		template<typename UnaryFunction>
		void call_for_each_boundary_edge(UnaryFunction f)
		{
			util::for_each_in_tuple(_boundary_edges.data,
				[f](auto& edge_vector) {
				std::for_each(edge_vector.begin(), edge_vector.end(),
					f);
			});
		}*/
	};
}

#endif
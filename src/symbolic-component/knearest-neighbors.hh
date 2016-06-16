// Copyright (c) 2014, LAAS-CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-manipulation.
// hpp-manipulation is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-manipulation is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-manipulation. If not, see <http://www.gnu.org/licenses/>.

#ifndef HPP_MANIPULATION_SYMBOLIC_COMPONENT_KNEAREST_NEIGHBORS_HH
#define HPP_MANIPULATION_SYMBOLIC_COMPONENT_KNEAREST_NEIGHBORS_HH

#include "hpp/manipulation/symbolic-component.hh"

#include <hpp/core/distance.hh>

namespace hpp {
  namespace manipulation {
    namespace symbolicComponent {
    class HPP_MANIPULATION_LOCAL KNearestNeighbors :
      public SymbolicComponentVisitor
    {
      private:
        typedef std::pair <value_type, RoadmapNodePtr_t> DistAndNode_t;
        struct DistAndNodeComp_t {
          bool operator () (const DistAndNode_t& r,
              const DistAndNode_t& l) {
            return r.first < l.first;
          }
        };
        typedef std::priority_queue <DistAndNode_t, std::vector <DistAndNode_t>,
                DistAndNodeComp_t > Queue_t;

      public:
        void visit (const SymbolicComponentPtr_t from)
        {
          while (!nearest_.empty()) nearest_.pop();

          const core::Distance& dist = *from->roadmap()->distance();
          const Configuration_t& q = *q_;
          const RoadmapNodes_t& nodes = from->nodes();
          for (RoadmapNodes_t::const_iterator _node = nodes.begin();
              _node != nodes.end (); ++_node) {
            const RoadmapNodePtr_t& node = *_node;
            value_type d = dist (*node->configuration (), q);
            if (nearest_.size () < k_)
              nearest_.push (DistAndNode_t (d, node));
            else if (nearest_.top().first > d) {
              nearest_.pop ();
              nearest_.push (DistAndNode_t (d, node));
            }
          }
        }

        /// \warning the content of the last visit is cleared.
        RoadmapNodes_t nodes ()
        {
          RoadmapNodes_t nodes(nearest_.size());
          std::size_t i = 0;
          while (!nearest_.empty()) {
            nodes[i] = nearest_.top().second;
            nearest_.pop ();
            ++i;
          }
          return nodes;
        }

        KNearestNeighbors (const std::size_t& k, const ConfigurationPtr_t& q)
          : k_ (k), q_ (q) {}
        
      private:
        const std::size_t k_;
        ConfigurationPtr_t q_;
        Queue_t nearest_;
    };
    } // namespace symbolicComponent
  } // namespace manipulation
} // namespace hpp

#endif // HPP_MANIPULATION_SYMBOLIC_COMPONENT_KNEAREST_NEIGHBORS_HH

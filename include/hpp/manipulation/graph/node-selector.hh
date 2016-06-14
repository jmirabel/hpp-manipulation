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

#ifndef HPP_MANIPULATION_GRAPH_NODE_SELECTOR_HH
# define HPP_MANIPULATION_GRAPH_NODE_SELECTOR_HH

#include "hpp/manipulation/config.hh"
#include "hpp/manipulation/fwd.hh"
#include "hpp/manipulation/graph/graph.hh"
#include "hpp/manipulation/graph/node.hh"

namespace hpp {
  namespace manipulation {
    namespace graph {
      /// This class is used to get the state of a configuration. States have to
      /// be ordered because a configuration can be in several states.
      class HPP_MANIPULATION_DLLAPI NodeSelector : public GraphComponent
      {
        public:
          /// Create a new NodeSelector.
          static NodeSelectorPtr_t create(const std::string& name);

          /// Create an empty node
          NodePtr_t createNode (const std::string& name, bool waypoint = false,
              const int w = 0);

          /// Returns the state of a configuration.
          NodePtr_t getNode(ConfigurationIn_t config) const;

          /// Returns the state of a roadmap node
          NodePtr_t getNode(RoadmapNodePtr_t node) const;

          /// Returns a list of all the nodes
          Nodes_t getNodes () const;

          /// Select randomly an outgoing edge of the given node.
          virtual EdgePtr_t chooseEdge(RoadmapNodePtr_t from) const;

          /// Should never be called.
          void addNumericalConstraint (
              const core::NumericalConstraintPtr_t& /* function */,
              const SizeIntervals_t& /* passiveDofs */ = SizeIntervals_t ())
          {
            HPP_THROW_EXCEPTION (Bad_function_call, "This component does not have constraints.");
          }

          /// Should never be called.
          void addLockedJointConstraint
	    (const core::LockedJoint& /* constraint */)
          {
            HPP_THROW_EXCEPTION (Bad_function_call, "This component does not have constraints.");
          }

          /// Print the object in a stream.
          virtual std::ostream& dotPrint (std::ostream& os, dot::DrawingAttributes da = dot::DrawingAttributes ()) const;

        protected:
          /// Initialization of the object.
          void init (const NodeSelectorPtr_t& weak);

          /// Constructor
          NodeSelector (const std::string& name) : GraphComponent (name)
          {}

          /// Print the object in a stream.
          virtual std::ostream& print (std::ostream& os) const;

          /// List of the states of one end-effector, ordered by priority.
          typedef std::pair <int, NodePtr_t> WeighedNode_t;
          typedef std::list <WeighedNode_t> WeighedNodes_t;
          WeighedNodes_t orderedStates_;
          WeighedNodes_t waypoints_;

        private:
          /// Weak pointer to itself.
          NodeSelectorPtr_t wkPtr_;
      }; // Class NodeSelector
    } // namespace graph
  } // namespace manipulation
} // namespace hpp

#endif // HPP_MANIPULATION_GRAPH_NODE_SELECTOR_HH

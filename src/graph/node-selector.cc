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

#include <hpp/core/node.hh>

#include <hpp/model/configuration.hh>
#include "hpp/manipulation/roadmap-node.hh"
#include "hpp/manipulation/graph/node-selector.hh"

#include <cstdlib>

namespace hpp {
  namespace manipulation {
    namespace graph {
      NodeSelectorPtr_t NodeSelector::create(const std::string& name)
      {
        NodeSelector* ptr = new NodeSelector (name);
        NodeSelectorPtr_t shPtr (ptr);
        ptr->init (shPtr);
        return shPtr;
      }

      void NodeSelector::init (const NodeSelectorPtr_t& weak)
      {
        GraphComponent::init (weak);
        wkPtr_ = weak;
      }

      NodePtr_t NodeSelector::createNode (const std::string& name,
          bool waypoint, const int w)
      {
        NodePtr_t newNode = Node::create (name);
        newNode->nodeSelector(wkPtr_);
        newNode->parentGraph(graph_);
        newNode->isWaypoint (waypoint);
        WeighedNodes_t& nodes_ = (waypoint?waypoints_:orderedStates_);
        bool found = false;
        for (WeighedNodes_t::iterator it = nodes_.begin();
            it != nodes_.end (); ++it) {
          if (it->first < w) {
            nodes_.insert (it, WeighedNode_t(w,newNode));
            found = true;
            break;
          }
        }
        if (!found) 
          nodes_.push_back (WeighedNode_t(w,newNode));
        return newNode;
      }

      Nodes_t NodeSelector::getNodes () const
      {
        Nodes_t ret;
        for (WeighedNodes_t::const_iterator it = orderedStates_.begin();
            it != orderedStates_.end (); ++it)
          ret.push_back (it->second);
        return ret;
      }

      NodePtr_t NodeSelector::getNode(ConfigurationIn_t config) const
      {
        for (WeighedNodes_t::const_iterator it = orderedStates_.begin();
	     orderedStates_.end() != it; ++it) {
          if (it->second->contains(config))
            return it->second;
	}
	std::stringstream oss;
	oss << "A configuration has no node:" << model::displayConfig (config);
	throw std::logic_error (oss.str ());
      }

      NodePtr_t NodeSelector::getNode(RoadmapNodePtr_t node) const
      {
        NodePtr_t n;
        switch (node->cachingSystem ()) {
          case RoadmapNode::CACHE_UP_TO_DATE:
            n = node->graphNode ();
            break;
          case RoadmapNode::CACHE_DISABLED:
          case RoadmapNode::CACHE_NEED_UPDATE:
            n = getNode (*(node->configuration ()));
            node->graphNode (n);
            break;
          default:
            n = getNode (*(node->configuration ()));
            hppDout (error, "Unimplemented caching system.");
            break;
        }
        return n;
      }

      EdgePtr_t NodeSelector::chooseEdge(RoadmapNodePtr_t from) const
      {
        NodePtr_t node = getNode (from);
        const Neighbors_t neighborPicker = node->neighbors();
        if (neighborPicker.totalWeight () == 0) {
          return EdgePtr_t ();
        }
        return neighborPicker ();
      }

      std::ostream& NodeSelector::dotPrint (std::ostream& os, dot::DrawingAttributes) const
      {
        for (WeighedNodes_t::const_iterator it = orderedStates_.begin();
            orderedStates_.end() != it; ++it)
          it->second->dotPrint (os);
        return os;
      }

      std::ostream& NodeSelector::print (std::ostream& os) const
      {
        os << "|-- ";
        GraphComponent::print (os) << std::endl;
        for (WeighedNodes_t::const_iterator it = orderedStates_.begin();
            orderedStates_.end() != it; ++it)
          os << it->first << " " << *it->second;
        return os;
      }
    } // namespace graph
  } // namespace manipulation
} // namespace hpp

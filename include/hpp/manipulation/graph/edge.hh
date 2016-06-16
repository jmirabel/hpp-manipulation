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

#ifndef HPP_MANIPULATION_GRAPH_EDGE_HH
# define HPP_MANIPULATION_GRAPH_EDGE_HH

#include <hpp/core/constraint-set.hh>
#include <hpp/core/steering-method.hh>
#include <hpp/core/relative-motion.hh>
#include <hpp/core/path.hh>

#include "hpp/manipulation/config.hh"
#include "hpp/manipulation/fwd.hh"
#include "hpp/manipulation/graph/graph.hh"

namespace hpp {
  namespace manipulation {
    namespace graph {
      /// \cond
      /// Cache mechanism that enable const-correctness of member functions.
      template <typename C>
        class HPP_MANIPULATION_LOCAL Cache
      {
        public:
          void set (const C& c)
          {
            c_ = c;
          }

          operator bool() const
          {
            return (bool)c_;
          }

          const C& get () const
          {
            return c_;
          }

        private:
          C c_;
      };
      /// \endcond

      /// \addtogroup constraint_graph
      /// \{

      /// Transition between two nodes of a constraint graph
      ///
      /// An edge stores two types of constraints.

      /// \li <b> Path constraints </b> should be safisfied by paths belonging
      /// to the edge. Along any path, the right hand side of the
      /// constraint is constant, but can differ between paths. For
      /// instance if an edge represents a transit path of a robot
      /// that can grasp an object, the right hand side of the
      /// constraint represents the position of the object. Along any
      /// transit path, the object does not move, but for different paths
      /// the object can be at different positions.
      /// \sa method pathConstraint.
      /// \li <b> Configuration constraints </b> are constraints that
      /// configurations in the destination node should satisfy and
      /// the constraints that paths should satisfy. For instance, if
      /// the edge links a node where the robot does not hold the
      /// object to a node where the robot holds the object, the
      /// configuration constraints represent a fixed relative
      /// position of the object with respect to the gripper and a
      /// stable position of the object. Configuration constraints are
      /// necessary to generate a configuration in the destination
      /// node of the edge that is reachable from a given
      /// configuration in the start node by an admissible path.
      class HPP_MANIPULATION_DLLAPI Edge : public GraphComponent
      {
        public:
          typedef core::RelativeMotion RelativeMotion;

          /// Destructor
          virtual ~Edge ();

          /// Create a new empty Edge.
          static EdgePtr_t create
	    (const std::string& name,
	     const GraphWkPtr_t& graph,
	     const NodeWkPtr_t& from,
	     const NodeWkPtr_t& to);

	  /// Apply edge constraint
	  ///
	  /// \param nnear node containing the configuration defining the right
	  ///        hand side of the edge constraint,
	  /// \param[in,out] q configuration to which the edge constraint is
	  ///                applied.
	  ///
	  /// \sa hpp::core::ConfigProjector::rightHandSideFromConfig
          virtual bool applyConstraints (core::NodePtr_t nnear, ConfigurationOut_t q) const;
	  /// Apply edge constraint
	  ///
	  /// \param qoffset configuration defining the right hand side of the
	  ///        edge constraint,
	  /// \param[in,out] q configuration to which the edge constraint is
	  ///                applied.
	  ///
	  /// \sa hpp::core::ConfigProjector::rightHandSideFromConfig
          virtual bool applyConstraints (ConfigurationIn_t qoffset, ConfigurationOut_t q) const;

          virtual bool canConnect (ConfigurationIn_t q1, ConfigurationIn_t q2) const;

          virtual bool build (core::PathPtr_t& path, ConfigurationIn_t q1,
              ConfigurationIn_t q2) const;

          /// Get the destination
          NodePtr_t to () const;

          /// Get the origin
          NodePtr_t from () const;

          /// Get the node in which path is.
          NodePtr_t node () const
          {
            return node_.lock();
          }

          void node (NodePtr_t node)
          {
            node_ = node;
          }

          /// \deprecated use node(NodePtr_t) instead.
          void isInNodeFrom (bool iinf) HPP_MANIPULATION_DEPRECATED
          {
            if (iinf) node_ = from_;
            else      node_ = to_;
          }

          /// \deprecated see NodePtr_t node() const
          bool isInNodeFrom () const HPP_MANIPULATION_DEPRECATED
          {
            return node_.lock() == from_.lock();
          }

	  /// Get steering method associated to the edge.
	  const core::SteeringMethodPtr_t& steeringMethod () const
	  {
	    return steeringMethod_->get();
	  }

	  /// Get path validation associated to the edge.
	  const core::PathValidationPtr_t& pathValidation () const
	  {
	    return pathValidation_->get();
	  }

          const RelativeMotion::matrix_type& relativeMotion () const
          {
            return relMotion_;
          }

          /// Get direction of the path compare to the edge
          /// \return true is reverse
          virtual bool direction (const core::PathPtr_t& path) const;

          /// Populate a ConfigProjector with constraints required to generate
          /// a path at the intersection of two edges.
          virtual bool intersectionConstraint (const EdgePtr_t& other,
              ConfigProjectorPtr_t projector) const;

          /// Print the object in a stream.
          virtual std::ostream& dotPrint (std::ostream& os, dot::DrawingAttributes da = dot::DrawingAttributes ()) const;

          /// Constraint of the destination node and of the path
          ConstraintSetPtr_t configConstraint() const;

          void setShort (bool isShort) {
            isShort_ = isShort;
          }

          bool isShort () const {
            return isShort_;
          }

        protected:
          /// Initialization of the object.
          void init (const EdgeWkPtr_t& weak, const GraphWkPtr_t& graph, const NodeWkPtr_t& from,
              const NodeWkPtr_t& to);

          /// Constructor
          Edge (const std::string& name);

          /// Constraint to project a path.
          /// \return The initialized constraint.
          ConstraintSetPtr_t pathConstraint() const;

          virtual ConstraintSetPtr_t buildConfigConstraint() const;

          virtual ConstraintSetPtr_t buildPathConstraint() const;

          /// Print the object in a stream.
          virtual std::ostream& print (std::ostream& os) const;

          bool isShort_;

        private:
          typedef Cache < ConstraintSetPtr_t > Constraint_t;
          typedef Cache < core::SteeringMethodPtr_t > SteeringMethod_t;
          typedef Cache < core::PathValidationPtr_t > PathValidation_t;

          /// See pathConstraint member function.
          Constraint_t* pathConstraints_;

          /// Constraint ensuring that a q_proj will be in to_ and in the
          /// same leaf of to_ as the configuration used for initialization.
          Constraint_t* configConstraints_;

          /// The two ends of the transition.
          NodeWkPtr_t from_, to_;

          /// True if this path is in node from, False if in node to
          NodeWkPtr_t node_;

	  /// Steering method used to create paths associated to the edge
	  SteeringMethod_t* steeringMethod_;

	  /// Path validation associated to the edge
          mutable RelativeMotion::matrix_type relMotion_;
	  PathValidation_t* pathValidation_;

          /// Weak pointer to itself.
          EdgeWkPtr_t wkPtr_;

          friend class Graph;
      }; // class Edge

      /// Edge with waypoint.
      /// Waypoints are handled recursively, i.e.\ class WaypointEdge contains only a
      /// Node and an Edge, the second Edge being itself.
      /// In this package, the Node in a WaypointEdge is semantically different from other Node
      /// because it does not correspond to a state with different manipulation rules. It has
      /// the same rules as another Node (either Edge::from() or Edge::to()).
      ///
      /// Semantically, a waypoint Node is fully part of the WaypointEdge. When a corresponding path
      /// reaches it, no planning is required to know what to do next. To the contrary, when a path reaches
      /// Edge::from() or Edge::to(), there may be several possibilities.
      ///
      /// \note
      ///   Implementation details: let's say, between the two nodes \f$N_f\f$ and \f$N_t\f$,
      ///   two waypoints are required:
      ///   \f$ N_f \xrightarrow{e_0} n_0 \xrightarrow{e_1} n_1 \xrightarrow{E} N_t\f$.
      ///   The outmost WaypointEdg contains:
      ///   \li from: \f$N_f\f$,
      ///   \li to: \f$N_t\f$,
      ///   \li constraints: those of edge \f$E\f$,
      ///   \li waypoint: \f$(E_1, n_1)\f$.
      ///
      ///   where \f$E_1\f$ is an instance of class WaypointEdge containing:
      ///   \li from: \f$N_f\f$,
      ///   \li to: \f$n_1\f$,
      ///   \li constraints: those of edge \f$e_1\f$,
      ///   \li waypoint: \f$(e_0, n_0)\f$.
      class HPP_MANIPULATION_DLLAPI WaypointEdge : public Edge
      {
        public:
          /// Create a new WaypointEdge.
	static WaypointEdgePtr_t create
	  (const std::string& name,
	   const GraphWkPtr_t& graph, const NodeWkPtr_t& from,
	   const NodeWkPtr_t& to);

          virtual bool direction (const core::PathPtr_t& path) const;

          virtual bool canConnect (ConfigurationIn_t q1, ConfigurationIn_t q2) const;

          virtual bool build (core::PathPtr_t& path, ConfigurationIn_t q1, ConfigurationIn_t q2) const;

          virtual bool applyConstraints (ConfigurationIn_t qoffset, ConfigurationOut_t q) const;

          /// Return the inner waypoint.
          /// \param EdgeType is either Edge or WaypointEdge
          template <class EdgeType>
          boost::shared_ptr <EdgeType> waypoint (const std::size_t index) const;

          /// Print the object in a stream.
          virtual std::ostream& dotPrint (std::ostream& os, dot::DrawingAttributes da = dot::DrawingAttributes ()) const;

          /// Set the number of waypoints
          void nbWaypoints (const size_type number);

          std::size_t nbWaypoints () const
          {
            return waypoints_.size ();
          }

          /// Set waypoint index with wEdge and wTo.
          /// \param wTo is the destination node of wEdge
          void setWaypoint (const std::size_t index, const EdgePtr_t wEdge, const NodePtr_t wTo);

          /// Get the node in which path after the waypoint is.
          NodePtr_t node () const;

        protected:
	  WaypointEdge (const std::string& name) :
	    Edge (name)
	    {
	    }
          /// Initialization of the object.
          void init (const WaypointEdgeWkPtr_t& weak, const GraphWkPtr_t& graph, const NodeWkPtr_t& from,
              const NodeWkPtr_t& to);

          /// Print the object in a stream.
          virtual std::ostream& print (std::ostream& os) const;

        private:
          typedef std::pair < EdgePtr_t, NodePtr_t > Waypoint_t;
          typedef std::vector <Waypoint_t> Waypoints_t;

          Waypoints_t waypoints_;

          mutable matrix_t configs_;
          mutable Configuration_t result_;

          WaypointEdgeWkPtr_t wkPtr_;
      }; // class WaypointEdge

      /// Edge that find intersection of level set.
      class HPP_MANIPULATION_DLLAPI LevelSetEdge : public Edge
      {
        public:
          virtual ~LevelSetEdge ();

          /// Create a new LevelSetEdge.
          static LevelSetEdgePtr_t create
	    (const std::string& name,
	     const GraphWkPtr_t& graph, const NodeWkPtr_t& from,
	     const NodeWkPtr_t& to);

          virtual bool applyConstraints (ConfigurationIn_t qoffset, ConfigurationOut_t q) const;

          virtual bool applyConstraints (core::NodePtr_t n_offset, ConfigurationOut_t q) const;

          virtual ConstraintSetPtr_t buildConfigConstraint() const;

          void buildHistogram ();

          LeafHistogramPtr_t histogram () const;

          /// \name Foliation definition
          /// \{

          /// Insert a NumericalConstraint that parametrizes the foliation
          void insertParamConstraint (const NumericalConstraintPtr_t& nm,
              const SizeIntervals_t& passiveDofs = SizeIntervals_t ());

          void insertParamConstraint (const DifferentiableFunctionPtr_t function, const ComparisonTypePtr_t ineq)
            HPP_MANIPULATION_DEPRECATED;

          /// Insert a LockedJoint that parametrizes the foliation
          void insertParamConstraint (const LockedJointPtr_t lockedJoint);

          /// Insert a NumericalConstraint that defines the foliation
          ///
          /// The manifold represented the foliation is defined by this
          /// constraints.
          void insertConditionConstraint (const NumericalConstraintPtr_t& nm,
              const SizeIntervals_t& passiveDofs = SizeIntervals_t ());

          /// Insert a LockedJoint that defines the foliation
          ///
          /// The manifold represented the foliation is defined by this
          /// constraints.
          void insertConditionConstraint (const LockedJointPtr_t lockedJoint);
          /// \}

          /// Print the object in a stream.
          virtual std::ostream& dotPrint (std::ostream& os, dot::DrawingAttributes da = dot::DrawingAttributes ()) const;

        protected:
          /// Initialization of the object.
          void init (const LevelSetEdgeWkPtr_t& weak, const GraphWkPtr_t& graph, const NodeWkPtr_t& from,
              const NodeWkPtr_t& to);

	  LevelSetEdge (const std::string& name);

          /// Print the object in a stream.
          virtual std::ostream& print (std::ostream& os) const;

          /// Populate DrawingAttributes tooltip
          virtual void populateTooltip (dot::Tooltip& tp) const;

        private:
          bool applyConstraintsWithOffset (ConfigurationIn_t qoffset,
              ConfigurationIn_t qlevelset, ConfigurationOut_t q) const;

          // Parametrizer
          // NumericalConstraints_t
          NumericalConstraints_t paramNumericalConstraints_;
          IntervalsContainer_t paramPassiveDofs_;
          // LockedJoints_t
          LockedJoints_t paramLockedJoints_;

          // Condition
          // NumericalConstraints_t
          NumericalConstraints_t condNumericalConstraints_;
          IntervalsContainer_t condPassiveDofs_;
          // LockedJoints_t
          LockedJoints_t condLockedJoints_;

          /// This histogram will be used to find a good level set.
          LeafHistogramPtr_t hist_;

          LevelSetEdgeWkPtr_t wkPtr_;
      }; // class LevelSetEdge

      /// \}
    } // namespace graph
  } // namespace manipulation
} // namespace hpp

#endif // HPP_MANIPULATION_GRAPH_EDGE_HH

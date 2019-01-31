// Copyright (c) 2019, Joseph Mirabel
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

#include <hpp/constraints/implicit.hh>
#include <hpp/constraints/differentiable-function.hh>

#include <hpp/core/path-vector.hh>
#include <hpp/core/plugin.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/problem-solver.hh>
#include <hpp/core/config-projector.hh>

#include <hpp/manipulation/steering-method/end-effector-trajectory.hh>

namespace hpp {
  namespace manipulation {
    using core::Parameter;
    using steeringMethod::EndEffectorTrajectory;
    using steeringMethod::EndEffectorTrajectoryPtr_t;

    class SetEndEffectorTrajectory : public core::Command
    {
      public:
        SetEndEffectorTrajectory ()
          : Command ("Setup EndEffectorTrajectory steering method.\n"
              "Arguments are:\n"
              "- constraint name\n"
              "- path index\n"
              "- start of time interval in path\n"
              "- end of time interval in path\n"
              "- True if path output space is SE3, false otherwise.\n"
              )
        {
          valueTypes (core::list_of
              (Parameter::STRING)
              (Parameter::INT)
              (Parameter::FLOAT)
              (Parameter::FLOAT)
              (Parameter::BOOL));
        }

      protected:
        Parameter doExecute (core::ProblemSolverPtr_t ps)
        {
          std::string cname = getParameterValues()[0].stringValue();
          size_type iPath   = getParameterValues()[1].intValue();
          value_type start  = getParameterValues()[2].floatValue();
          value_type end    = getParameterValues()[3].floatValue();
          bool se3output    = getParameterValues()[4].boolValue();

#if 1
          constraints::ImplicitPtr_t constraint = ps->numericalConstraint (cname);
          if (!constraint)
            throw std::invalid_argument ("Constraint not found.");
#else
          core::ConstraintSetPtr_t c (ps->constraints());
          if (!c || !c->configProjector())
            throw std::logic_error ("Problem constraint must have a ConfigProjector");
          ConfigProjectorPtr_t cp ();

          constraints::ImplicitPtr_t constraint;
          const core::NumericalConstraints_t& ncs = c->configProjector()->numericalConstraints();
          for (std::size_t i = 0; i < ncs.size(); ++i) {
            if (ncs[i]->function().name() == cname) {
              constraint = ncs[i];
              break;
            }
          }
#endif

          if (iPath < 0 || (std::size_t)iPath >= ps->paths().size())
            throw std::invalid_argument ("Path index out of range");
          core::PathPtr_t path = ps->paths()[iPath];
          core::interval_t tr = path->timeRange ();
          if (tr.first != start || tr.second != end) {
            tr = core::interval_t(start,end);
            path = path->extract (tr);
          }

          EndEffectorTrajectoryPtr_t sm = HPP_DYNAMIC_PTR_CAST (EndEffectorTrajectory,
              ps->problem()->steeringMethod ());
          if (!sm)
            throw std::invalid_argument ("Steering method is not of type EndEffectorTrajectory");

          sm->trajectoryConstraint (constraint->copy());
          sm->trajectory (path, se3output);

          return Parameter();
        }
    };

    class EndEffectorTrajectoryPlugin : public core::ProblemSolverPlugin
    {
      public:
        EndEffectorTrajectoryPlugin ()
          : ProblemSolverPlugin ("EndEffectorTrajectoryPlugin", "0.0")
        {}

      protected:
        virtual bool impl_initialize (core::ProblemSolverPtr_t ps)
        {
          ps->steeringMethods.add ("EndEffectorTrajectory", EndEffectorTrajectory::create);
          ps->commands.add ("setupEndEffectorTrajectory",
              core::CommandPtr_t (new SetEndEffectorTrajectory()));

          return true;
        }
    };
  } // namespace manipulation
} // namespace hpp

HPP_CORE_DEFINE_PLUGIN(hpp::manipulation::EndEffectorTrajectoryPlugin)

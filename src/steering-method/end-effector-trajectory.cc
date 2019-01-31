// Copyright (c) 2019, LAAS-CNRS
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

#include <hpp/manipulation/steering-method/end-effector-trajectory.hh>

#include <hpp/util/indent.hh>

#include <hpp/pinocchio/util.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/urdf/util.hh>

#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/implicit.hh>

#include <hpp/core/config-projector.hh>
#include <hpp/core/path.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/straight-path.hh>

namespace hpp {
  namespace manipulation {
    namespace steeringMethod {
      namespace {
        /// Class only for testing purpose.
        class Interpolation : public constraints::DifferentiableFunction
        {
          public:
            Interpolation (const matrix_t& points,
                const vector_t& times)
              : DifferentiableFunction (1, 1, points.cols()),
              points_ (points),
              times_ (times)
            {}

          protected:
            void impl_compute (core::LiegroupElementRef result, vectorIn_t arg) const
            {
              value_type t = arg[0];
              size_type i = rank (t);
              value_type u = (t - times_[i]) / (times_[i+1] - times_[i]);
              result.vector() = (1-u) * points_.row(i) + u * points_.row(i+1);
            }

            void impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const
            {
              value_type t = arg[0];
              size_type i = rank (t);
              jacobian = points_.row(i+1) - points_.row(i);
            }

            std::ostream& print (std::ostream& os) const
            {
              return os << "Interpolation: " << one_line (times_);
            }

          private:
            size_type rank (value_type t) const
            {
              for (size_type i = 1; i < times_.size(); ++i)
                if (t < times_[i]) return i-1; // between i-1 and i
              return times_.size() - 2; // between i-2 and i-1
            }

            matrix_t points_;
            vector_t times_;
        };

        template <bool SE3>
        class FunctionFromPath : public constraints::DifferentiableFunction
        {
          public:
            FunctionFromPath (const PathPtr_t& p)
              : DifferentiableFunction (1, 1, SE3 ? LiegroupSpace::SE3() : LiegroupSpace::Rn(p->outputSize())),
              path_ (p)
            {
              assert (!SE3 ||
                  (p->outputSize() == 7 && p->outputDerivativeSize() == 6));
            }

            const PathPtr_t& path () const
            {
              return path_;
            }

            std::ostream& print (std::ostream& os) const
            {
              return os
                << (SE3 ? "FunctionFromSE3Path: " : "FunctionFromPath: ")
                << path_->timeRange().first << ", "
                << path_->timeRange().second
                ;
            }

          protected:
            void impl_compute (core::LiegroupElementRef result, vectorIn_t arg) const
            {
              bool success = (*path_) (result.vector(), arg[0]);
              if (!success) {
                hppDout (warning, "Failed to evaluate path at param " << arg[0]
                    << incindent << iendl << *path_ << decindent);
              }
            }

            void impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const
            {
              path_->derivative (jacobian.col(0), arg[0], 1);
            }

          private:
            PathPtr_t path_;
        };
      }

      EndEffectorTrajectoryPtr_t EndEffectorTrajectory::createWithGuess
        (const core::Problem& problem)
      {
        matrix_t points = problem.getParameter ("EndEffectorTrajectory/points").matrixValue();
        vector_t times = problem.getParameter ("EndEffectorTrajectory/times").vectorValue();
        std::string cname = problem.getParameter ("EndEffectorTrajectory/constraint").stringValue();

        core::ConstraintSetPtr_t c (problem.constraints());
        if (!c || !c->configProjector())
          throw std::logic_error ("Problem constraint must have a ConfigProjector");
        ConfigProjectorPtr_t cp (c->configProjector());

        ImplicitPtr_t ic;
        const core::NumericalConstraints_t& ncs = cp->numericalConstraints();
        for (std::size_t i = 0; i < ncs.size(); ++i) {
          if (ncs[i]->function().name() == cname) {
            ic = ncs[i];
            break;
          }
        }
        if (!ic)
          throw std::logic_error ("Constraint not found");
        if (times[0] != 0)
          throw std::logic_error ("Times[0] must be 0");
        for (size_type i = 1; i < times.size(); ++i)
          if (times[i-1] >= times[i])
            throw std::logic_error ("Times must be strictly increasing");
        if (points.rows() != times.size())
          throw std::logic_error ("Points and times size mismatch.");
        if (points.cols() != ic->rhsSize())
          throw std::logic_error ("Points size does not match constraint right hand side.");

        DifferentiableFunctionPtr_t path (new Interpolation (points, times));

        EndEffectorTrajectoryPtr_t sm = create (problem);
        sm->trajectory (path, interval_t(times[0], times[times.size()-1]));
        sm->trajectoryConstraint (ic);

        return sm;
      }

      void EndEffectorTrajectory::trajectoryConstraint (const constraints::ImplicitPtr_t& ic)
      {
        constraint_ = ic;
        if (eeTraj_) trajectory (eeTraj_, timeRange_);
      }

      void EndEffectorTrajectory::trajectory (const PathPtr_t& eeTraj, bool se3Output)
      {
        if (se3Output)
          trajectory (DifferentiableFunctionPtr_t
              (new FunctionFromPath <true > (eeTraj)), eeTraj->timeRange());
        else
          trajectory (DifferentiableFunctionPtr_t
              (new FunctionFromPath <false> (eeTraj)), eeTraj->timeRange());
      }

      void EndEffectorTrajectory::trajectory (const DifferentiableFunctionPtr_t& eeTraj, const interval_t& tr)
      {
        assert (eeTraj->inputSize() == 1);
        assert (eeTraj->inputDerivativeSize() == 1);
        eeTraj_ = eeTraj;
        timeRange_ = tr;

        if (!constraint_) return;

        constraint_->rightHandSideFunction (eeTraj_);
      }

      PathPtr_t EndEffectorTrajectory::impl_compute (ConfigurationIn_t q1,
          ConfigurationIn_t q2) const
      {
        try {
        if (!eeTraj_) throw std::logic_error ("EndEffectorTrajectory not initialized.");
        // Update (or check) the constraints
        core::ConstraintSetPtr_t c (constraints());
        if (!c || !c->configProjector()) {
          throw std::logic_error ("EndEffectorTrajectory steering method must "
              "have a ConfigProjector");
        }
        ConfigProjectorPtr_t cp (c->configProjector());

        const core::NumericalConstraints_t& ncs = cp->numericalConstraints();
        bool ok = false;
        for (std::size_t i = 0; i < ncs.size(); ++i) {
          if (ncs[i] == constraint_) {
            ok = true;
            break; // Same pointer
          }
          // Here, we do not check the right hand side on purpose.
          // if (*ncs[i] == *constraint_) {
          if (ncs[i]->functionPtr() == constraint_->functionPtr()
              && ncs[i]->comparisonType() == constraint_->comparisonType()) {
            ok = true;
            // TODO We should only modify the path constraint.
            // However, only the pointers to implicit constraints are copied
            // while we would like the implicit constraints to be copied as well.
            ncs[i]->rightHandSideFunction (eeTraj_);
            break; // logically identical
          }
        }
        if (!ok) {
          HPP_THROW (std::logic_error, "EndEffectorTrajectory could not find "
              "constraint " << constraint_->function());
        }

        return core::StraightPath::create (problem().robot(), q1, q2, timeRange_, c);
        } catch (const std::exception& e) {
          std::cout << timeRange_.first << ", " << timeRange_.second << '\n';
          if (eeTraj_)
            std::cout << (*eeTraj_) (vector_t::Constant(1,timeRange_.first )) << '\n'
                      << (*eeTraj_) (vector_t::Constant(1,timeRange_.second)) << std::endl;
          std::cout << *constraints() << std::endl;
          std::cout << e.what() << std::endl;
          throw;
        }
      }

      using core::Parameter;
      using core::ParameterDescription;

      HPP_START_PARAMETER_DECLARATION(EndEffectorTrajectory)
      core::Problem::declareParameter(ParameterDescription(Parameter::STRING,
            "EndEffectorTrajectory/constraint",
            "Name of the constraint in the problem constraint which represent"
            ));
      core::Problem::declareParameter(ParameterDescription(Parameter::MATRIX,
            "EndEffectorTrajectory/points",
            ""));
      core::Problem::declareParameter(ParameterDescription(Parameter::VECTOR,
            "EndEffectorTrajectory/times",
            ""));
      HPP_END_PARAMETER_DECLARATION(EndEffectorTrajectory)
    } // namespace steeringMethod
  } // namespace manipulation
} // namespace hpp

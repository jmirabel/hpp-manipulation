// Copyright (c) 2017, Joseph Mirabel
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

#include <hpp/constraints/differentiable-function.hh>

#include <hpp/manipulation/config.hh>

namespace hpp {
  namespace manipulation {
    namespace steeringMethod {
      typedef core::segment_t segment_t;
      namespace {
        std::string toStr (const segment_t& s)
        {
          std::ostringstream os;
          os << "[ " << s.first << ", " << s.first + s.second << " ]";
          return os.str();
        }
      }

      /// Apply the constraint on a subspace of the input space.
      /// i.e.: \f$ f (q_0, ... , q_n) = f_{inner} (q_k) \f$
      class HPP_MANIPULATION_LOCAL StateFunction :
        public constraints::DifferentiableFunction
      {
        public:
          typedef boost::shared_ptr<StateFunction> Ptr_t;
          StateFunction (const DifferentiableFunctionPtr_t& inner,
              const size_type& nArgs, const size_type& nDers,
              const segment_t& inArgs, const segment_t& inDers) :
            DifferentiableFunction (nArgs, nDers,
                inner->outputSpace(), inner->name() + " | " + toStr(inArgs)),
            inner_ (inner), sa_ (inArgs), sd_ (inDers)
          {
            activeParameters_.setConstant(false);
            activeParameters_.segment(sa_.first, sa_.second)
              = inner_->activeParameters();

            activeDerivativeParameters_.setConstant(false);
            activeDerivativeParameters_.segment(sd_.first, sd_.second)
              = inner_->activeDerivativeParameters();

            hppDout (info, inner_->name() << ": "
              << inner_->activeParameters().transpose());
          }

        protected:
          void impl_compute (LiegroupElement& y, vectorIn_t arg) const
          {
            inner_->value(y, arg.segment (sa_.first, sa_.second));
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t arg) const
          {
            inner_->jacobian(J.middleCols (sd_.first, sd_.second),
                arg.segment (sa_.first, sa_.second));
          }

          std::ostream& print (std::ostream& os) const
          {
            constraints::DifferentiableFunction::print(os);
            return os << incindent << iendl << *inner_ << decindent;
          }

          DifferentiableFunctionPtr_t inner_;
          const segment_t sa_, sd_;
      }; // class Function

      /// \f$ q_{out} = q_{in} \f$
      /// \todo Make this derive from constraints::AffineFunction
      class HPP_MANIPULATION_LOCAL Identity :
        public constraints::DifferentiableFunction
      {
        public:
          typedef boost::shared_ptr<Identity> Ptr_t;

          Identity (const LiegroupSpacePtr_t space, const std::string& name) :
            DifferentiableFunction (space->nq(), space->nv(), space, name) {}

        protected:
          void impl_compute (LiegroupElement& y, vectorIn_t arg) const
          {
            y.vector() = arg;
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t) const
          {
            J.setIdentity();
          }
      }; // class Function

      /// Compute the difference between the value of the function in two points.
      /// i.e.: \f$ f (q_0, ... , q_n) = f_{inner} (q_{left}) - f_{inner} (q_{right}) \f$
      class HPP_MANIPULATION_LOCAL EdgeFunction :
        public constraints::DifferentiableFunction
      {
        public:
          typedef boost::shared_ptr<EdgeFunction> Ptr_t;
          EdgeFunction (const DifferentiableFunctionPtr_t& inner,
              const size_type& nArgs, const size_type& nDers,
              const segment_t& lInArgs, const segment_t& lInDers,
              const segment_t& rInArgs, const segment_t& rInDers) :
            DifferentiableFunction (nArgs, nDers,
                LiegroupSpace::Rn(inner->outputSpace()->nv()),
                inner->name() + " | " + toStr(lInArgs) + " - " + toStr(rInArgs)),
            inner_ (inner),
            lsa_ (lInArgs), lsd_ (lInDers),
            rsa_ (rInArgs), rsd_ (rInDers),
            l_ (inner->outputSpace()), r_ (inner->outputSpace())
          {
            activeParameters_.setConstant(false);
            activeParameters_.segment(lsa_.first, lsa_.second)
              = inner_->activeParameters();
            activeParameters_.segment(rsa_.first, rsa_.second)
              = inner_->activeParameters();

            activeDerivativeParameters_.setConstant(false);
            activeDerivativeParameters_.segment(lsd_.first, lsd_.second)
              = inner_->activeDerivativeParameters();
            activeDerivativeParameters_.segment(rsd_.first, rsd_.second)
              = inner_->activeDerivativeParameters();
          }

        protected:
          void impl_compute (LiegroupElement& y, vectorIn_t arg) const
          {
            inner_->value(l_, arg.segment (lsa_.first, lsa_.second));
            inner_->value(r_, arg.segment (rsa_.first, rsa_.second));
            y.vector() = l_ - r_;
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t arg) const
          {
            inner_->jacobian(
                J.middleCols (lsd_.first, lsd_.second),
                arg.segment (lsa_.first, lsa_.second));
            inner_->jacobian(
                J.middleCols (rsd_.first, rsd_.second),
                arg.segment (rsa_.first, rsa_.second));
            J.middleCols (rsd_.first, rsd_.second) *= -1;
          }

          std::ostream& print (std::ostream& os) const
          {
            constraints::DifferentiableFunction::print(os);
            return os << incindent << iendl << *inner_ << decindent;
          }

          DifferentiableFunctionPtr_t inner_;
          const segment_t lsa_, lsd_;
          const segment_t rsa_, rsd_;

          mutable LiegroupElement l_, r_;
      }; // class Function

      /// Cost function: distances between the configurations.
      /// Weights are added for transition which should be short.
      class HPP_MANIPULATION_LOCAL CostFunction :
        public constraints::DifferentiableFunction
      {
        public:
          typedef boost::shared_ptr<CostFunction> Ptr_t;
          CostFunction (const core::DevicePtr_t& d, const size_type& N):
            DifferentiableFunction (N * d->configSize(), N * d->numberDof(),
                LiegroupSpace::R1(), "CostFunction"),
            robot_ (d),
            N_ (N),
            wdofs_ (vector_t::Ones(d->numberDof())),
            wcfgs_ (vector_t::Ones(N-1))
          {
            v[0] = vector_t(d->numberDof());
            v[1] = vector_t(d->numberDof());
          }

        protected:
          void impl_compute (LiegroupElement& y, vectorIn_t arg) const
          {
            size_type nq = robot_->configSize(),
                      i = nq,
                      j = 0;
            y.vector()[0] = 0;
            for (size_type k = 0; k < N_-1; ++k) {
              // v = q_{k+1} - q_{k}
              pinocchio::difference (robot_, arg.segment (j, nq), arg.segment (i, nq), v[0]);
              y.vector()[0] += wcfgs_[k] * v[0].cwiseAbs2().cwiseProduct(wdofs_).sum();
              j = i;
              i += nq;
            }
            y.vector()[0] /= 2;
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t arg) const
          {
            size_type nv = robot_->numberDof(),
                      nq = robot_->configSize(),
                      l = 0,
                      i = nq,
                      j = 0;

            int a = 0, b = 1;
            value_type wa, wb;
            for (size_type k = 0; k < N_; ++k) {
              a = (int) k%2;
              b = (a+1)%2;
              // v_a = q_{k+1} - q_{k}
              // v_b = q_{k} - q_{k-1}
              if (k == 0) {
                v[b].setZero();
                wa = wcfgs_[0];
                wb = 0;
              } else if (k == N_ - 1) {
                v[a].setZero();
                wa = 1;
                wb = wcfgs_[N_-2];
              } else {
                pinocchio::difference (robot_, arg.segment (j, nq), arg.segment (i, nq), v[a]);
                wa = wcfgs_[k  ];
                wb = wcfgs_[k-1];
              }

              J.middleCols (l, nv).noalias() = (wb * v[b] - wa * v[a]).cwiseProduct(wdofs_).transpose();

              j = i;
              i += nq;
              l += nv;
            }
          }

          core::DevicePtr_t robot_;
          const size_type N_;
          vector_t wdofs_,wcfgs_;

          mutable vector_t v[2];
      }; // class CostFunction
    } // namespace steeringMethod
  } // namespace manipulation
} // namespace hpp

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

#ifndef HPP_MANIPULATION_SYMBOLIC_COMPONENT_PRINT_HH
#define HPP_MANIPULATION_SYMBOLIC_COMPONENT_PRINT_HH

#include "hpp/manipulation/symbolic-component.hh"

namespace hpp {
  namespace manipulation {
    namespace symbolicComponent {
    class HPP_MANIPULATION_LOCAL Print :
      public SymbolicComponentVisitor
    {
      public:
        void visit (const SymbolicComponentPtr_t from) {
          out_ << "SymbolicComponent";
          generic(from);
          out_ << "\n";
        }
        void visit (const WeighedSymbolicComponentPtr_t from) {
          out_ << "WeighedSymbolicComponent";
          generic(from);
          out_ << ", P=(" << from->p_.transpose() << "), W=" << from->weight_;
        }

        Print (std::ostream& o) : out_ (o) {}
        
      private:
        void generic(const SymbolicComponentPtr_t& from) {
          out_ << "(" << from.get() << "), state=" << from->state()->name()
            << ", N=" << from->nodes().size();
        }
        std::ostream& out_;
    };
    } // namespace symbolicComponent
  } // namespace manipulation
} // namespace hpp

#endif // HPP_MANIPULATION_SYMBOLIC_COMPONENT_PRINT_HH

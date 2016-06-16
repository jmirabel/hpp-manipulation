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

#ifndef HPP_MANIPULATION_SYMBOLIC_COMPONENT_BFS_HH
#define HPP_MANIPULATION_SYMBOLIC_COMPONENT_BFS_HH

#include "hpp/manipulation/symbolic-component.hh"

#include <queue>

#include "print.hh"

namespace hpp {
  namespace manipulation {
    namespace symbolicComponent {
      class HPP_MANIPULATION_LOCAL Bfs :
        public SymbolicComponentVisitor
      {
        public:
          typedef std::list<SymbolicComponentPtr_t> path_type; 

          void visit (const SymbolicComponentPtr_t from) {
            from_ = from;
            parent_.clear();
            if (from == target_) {
              found_ = true;
              return;
            }

            std::queue<SymbolicComponentPtr_t> queue;
            queue.push(from);
            SymbolicComponentPtr_t current;
            found_ = false;
            std::vector<SymbolicComponentPtr_t> done;
            while (!queue.empty ()) {
              current = queue.front();
              queue.pop();
              assert (current == from_ || parent_.find(current) != parent_.end());
              const SymbolicComponents_t& scs = current->to();
              for (SymbolicComponents_t::const_iterator _sc = scs.begin();
                  _sc != scs.end(); ++_sc) {
                const SymbolicComponentPtr_t& next = *_sc;
                if (next == from_ || parent_.find(next) != parent_.end())
                  continue; // Already visited
                parent_ [next] = current;
                if (next == target_) {
                  found_ = true;
                  return;
                }
                queue.push (next);
              }
              done.push_back(current);
            }
          }

          path_type path() const
          {
            path_type path;
            if (!found_) return path;
            path.push_front(target_);
            if (from_ == target_) return path;

            std::map<SymbolicComponentPtr_t,SymbolicComponentPtr_t>::const_iterator
              _current;
            for (_current = parent_.find(target_);
                _current != parent_.end() && _current->second != from_;
                _current = parent_.find(_current->second))
              path.push_front(_current->second);
            if (_current == parent_.end()) return path_type();
            path.push_front(from_);
            return path;
          }

          Bfs (const SymbolicComponentPtr_t target)
            : target_ (target), found_(false) {}

          std::ostream& printPath (std::ostream& out) const
          {
            if (!found_) {
              out << "Path not found.";
              return out;
            }
            Print printer(out);
            path_type p = path();
            out << "===== Begin symbolic component path =====\n";
            for (path_type::const_iterator _sc = p.begin(); _sc != p.end(); ++_sc) {
              (*_sc)->accept(printer);
              out << "\n";
            }
            out << "===== End symbolic component path =====";
            return out;
          }

        private:
          SymbolicComponentPtr_t target_, from_;
          bool found_;
          std::map <SymbolicComponentPtr_t,SymbolicComponentPtr_t> parent_;
      };
    } // namespace symbolicComponent
  } // namespace manipulation
} // namespace hpp

#endif // HPP_MANIPULATION_SYMBOLIC_COMPONENT_BFS_HH

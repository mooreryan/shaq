# Copyright 2016 Ryan Moore
# Contact: moorer@udel.edu
#
# This file is part of Shaq.
#
# Shaq is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shaq is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Shaq.  If not, see <http://www.gnu.org/licenses/>.

class Read
  def self.k_plus_one_mers read, ksize
    k = ksize+1
    (0..read.length-k).each_with_index do |start, pos|
      yield read[start, k], pos
    end
  end
end

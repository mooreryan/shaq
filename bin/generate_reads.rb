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

readsize = 25
numreads = 100
dna = "i_love_to_eat_apple_pie_and_love_to_eat_other_fun_things"

start_range = (0..dna.length-readsize).to_a
reads = numreads.times.each_with_index do |n|
  start = start_range.sample

  puts [">read_#{n}", dna[start, readsize]].join "\n"
end

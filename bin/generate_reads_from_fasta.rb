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

require "parse_fasta"

readsize = 100
reads_per_seq = 3000
err_rate = 0.02

FastaFile.open(ARGV.first).each_record_fast do |head, seq|
  head.upcase!
  start_range = (0..seq.length-readsize).to_a
  reads = reads_per_seq.times.each_with_index do |n|
    start = start_range.sample

    new_read = ""

    seq[start, readsize].chars.each_with_index do |char, idx|
      if rand < err_rate
        new_char = (%w[a c t g] - %w[char]).sample
        new_read << new_char
      else
        new_read << char
      end
    end

    puts ">seq_#{head}_read_#{n}\n#{new_read}"
  end
end

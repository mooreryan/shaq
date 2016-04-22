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

require_relative "read"

class Graph
  attr_accessor :graph, :kmer_to_read, :incounts

  def initialize ksize
    @ksize = ksize
    @graph = {}
    @kmer_to_read = {}
    @incounts = Hash.new 0
  end

  def add_to_graph header, read
    k1, k2 = "", ""
    pos = -1
    Read.k_plus_one_mers(read, @ksize) do |kmer, start_pos|
      k1 = kmer[0..-2]
      k2 = kmer[1..-1]

      if @graph.has_key? k1
        @graph[k1] << k2
        @kmer_to_read[k1] << [header, start_pos]
      else
        @graph[k1] = Set.new [k2]
        @kmer_to_read[k1] = [[header, start_pos]]
      end

      @incounts[k2] += 1

      pos = start_pos
    end

    # make sure to add the last k-1mer
    if @graph.has_key? k2
      @kmer_to_read[k2] << [header, pos+1]
    else
      @kmer_to_read[k2] = [[header, pos+1]]
    end
  end

  def kmer_abundance
    @kmer_to_read.map.with_object(Hash.new) do |(kmer, posns), hash|
      hash[kmer] = posns.count
    end
  end
end

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

require "set"
require "fileutils"
require_relative "read"

class Graph
  attr_accessor :graph, :kmer_to_read, :contig_connections, :contigs

  def initialize ksize
    @ksize = ksize
    @graph = {}
    @kmer_to_read = {}
    @branch_point_kmers = Set.new

    #{ kmer => [ [ctg1 with kmer, ctg2 with kmer], ... ]
    @contig_connections = {}
    @contigs = []
  end

  def add header, read
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

  def make_contigs
    coming_from_branch_point = false
    kmer = String.new @graph.first.first

    next_kmer = ""
    prev_kmer = ""
    walk = String.new kmer

    c_num = 0

    while !next_kmer.nil?
      if @branch_point_kmers.include?(prev_kmer)
        if @contig_connections.has_key? prev_kmer
          @contig_connections[prev_kmer] << c_num
        else
          @contig_connections[prev_kmer] = Set.new [c_num]
        end
      end

      if @graph.has_key?(kmer) && !@graph[kmer].empty?
        if branch_point?(@graph[kmer]) || @branch_point_kmers.include?(kmer)
          @branch_point_kmers << kmer

          if @contig_connections.has_key? kmer
            @contig_connections[kmer] << c_num
          else
            @contig_connections[kmer] = Set.new [c_num]
          end

          @contigs[c_num] = walk
          c_num += 1

          next_kmer = String.new @graph[kmer].first
          @graph[kmer].delete next_kmer
          walk = String.new next_kmer
        else
          next_kmer = String.new @graph[kmer].first
          walk << next_kmer[-1]

          @graph[kmer].delete next_kmer
        end

        @graph.delete(kmer) if @graph[kmer].empty?
      elsif !@graph.has_key?(kmer) && !walk_is_done? # end of contig more to go
        @contigs[c_num] = walk
        c_num += 1

        next_kmer = String.new @graph.first.first
        walk = String.new next_kmer
      elsif !@graph.has_key?(kmer) && walk_is_done?
        @contigs[c_num] = walk
        c_num += 1

        break
      end

      prev_kmer = kmer
      kmer = next_kmer
    end
  end

  # @note Makes the directory if needed
  def print_fasta fname
    FileUtils.mkdir_p File.dirname(fname)

    File.open(fname, "w") do |f|
      @contigs.each_with_index do |contig, idx|
        f.puts ">contig_#{idx+1}\n#{contig}"
      end
    end
  end

  # @note Makes the directory if needed
  def print_contig_connections fname
    FileUtils.mkdir_p File.dirname(fname)

    File.open(fname, "w") do |f|
      @contig_connections.each do |kmer, contigs|
        contigs.to_a.combination(2).each do |(ctg1, ctg2)|
          f.puts ["contig_#{ctg1}", "contig_#{ctg2}", kmer].join "\t"
        end
      end
    end
  end

  private

  def walk_is_done?
    @graph.empty?
  end

  def branch_point? set
    set.count > 1
  end
end

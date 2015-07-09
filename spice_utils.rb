# $URL: https://repos.ser.asu.edu/svn/lroc/ruby_spice/tags/release_20140107.1/spice_utils.rb $
# # Copyright (C) 2013 by Arizona State University and Mark Robinson.
# All rights reserved.
#
# This file has been released under the modified BSD license.
# See COPYING in the distribution package for details.
#
# Author: Nick Estes <nme@ser.asu.edu>
#
# $Author: joe $
#   $Date: 2014-01-07 14:34:02 -0700 (Tue, 07 Jan 2014) $
#    $Rev: 13467 $
#
#This file contains pieces of the ruby spice library which are written in ruby.


#This class is used to parse out mk files
#LROC mk files as a whole wayyy blow past the 1000 file handle limit
#This pares it down to only the kernels which we need.
class SpiceKernelHelper
  attr_reader :base_kernels

  def initialize
    @base_kernels = []
    @detail_kernels = []
  end

  def parse file
    mkdata = open(file).read.gsub(/\s+/, " ").strip

    mkdata =~ /PATH_VALUES\s*=\s*\(\s*([^)]+?)\s*\)/
      path_values = $1.split(/\s*,\s*/).collect {|v| v.gsub "'", ''}

    mkdata =~ /PATH_SYMBOLS\s*=\s*\(\s*([^)]+?)\s*\)/
      path_symbols = $1.split(/\s*,\s*/).collect {|v| v.gsub "'", ''}

    mkdata =~ /KERNELS_TO_LOAD\s*=\s*\(\s*([^)]+?)\s*\)/
      raw_kernels = $1

    path_symbols.zip(path_values).each do |key, value|
      raw_kernels.gsub!("$#{key}", value)
    end

    kernels = raw_kernels.split(/\s*,\s*/).collect {|v| v.gsub "'", ''}.collect

    kernels.each do |kernel|
      if kernel =~ /\w\w\w\w\w\w?_(\d\d\d\d\d\d\d)_(\d\d\d\d\d\d\d)_.*?.(bc|bsp)$/
        start_doy = $1
        stop_doy = $2

        if start_doy == stop_doy
          @base_kernels << kernel
        else
          start_time = Time.gm(start_doy[0..3].to_i, 1, 1, 0, 0, 0) + (start_doy[4..6].to_i-1)*24*60*60
          stop_time = Time.gm(stop_doy[0..3].to_i, 1, 1, 0, 0, 0) + (stop_doy[4..6].to_i-1)*24*60*60

          @detail_kernels << [start_time..stop_time, kernel]
        end
      else
        @base_kernels << kernel
      end
    end
  end

  def detail_kernels time
    magic_sym = RUBY_VERSION =~ /^1\.8\./ ? :include? : :cover?
    @detail_kernels.find_all do |range, kernel|
      range.send(magic_sym, time) or range.send(magic_sym, time + 86400) or range.send(magic_sym, time - 86400)
    end.collect { |range, kernel| kernel }
  end
end

#We need to add some functionality to Range for wndifd
class Range
  def intersection(other)
    my_min, my_max = first, exclude_end? ? max : last
    other_min, other_max = other.first, other.exclude_end? ? other.max : other.last

    new_min = self === other_min ? other_min : other === my_min ? my_min : nil
    new_max = self === other_max ? other_max : other === my_max ? my_max : nil

    new_min && new_max ? new_min..new_max : nil
  end
end

#####
#
#Special Ruby implementations of select cspice functions
#
######
module Spice

  #wndifd(window_one, window_two)
  #Takes the difference between two 'ruby spice windows'
  #cspice functions which return spice windows return a 2d array of windows in ruby spice
  #that is: [[beg, end], ... , [beg, end]] 
  #This function will return everything that is in window_one, but NOT in window two
  def wndifd one,two
    unless (one.kind_of?(Array) and two.kind_of?(Array))
      raise ArgumentError
    end
    return one if two.length == 0
    return []  if one.length == 0

    ret = []
    result = []
    one.each do |a|
      a_rng = (a[0]..a[1]) #first
      r =  a_rng.intersection(0..two[0][0])
      ret << r unless r.nil?
      (0...two.size-1).each do |b|

        r = a_rng.intersection(two[b][1]..two[b+1][0])
        ret << r unless r.nil?

      end
      r = a_rng.intersection(two.last.last..a[1]) #last
      ret << r unless r.nil?
    end
    ret.each {|range| result << [range.begin, range.end]}
    result
  end

end

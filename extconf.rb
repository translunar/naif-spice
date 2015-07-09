#!/usr/bin/ruby
#
# $URL: https://repos.ser.asu.edu/svn/lroc/ruby_spice/tags/release_20140107.1/extconf.rb $
# 
# Copyright (C) 2013 by Arizona State University and Mark Robinson.
# All rights reserved.
#
# This file has been released under the modified BSD license.
# See COPYING in the distribution package for details.
#
# Author: Nick Estes <nme@ser.asu.edu>
#
# $Author: alicht $
#   $Date: 2013-07-31 16:30:44 -0700 (Wed, 31 Jul 2013) $
#    $Rev: 12001 $
#

# cspice can be compiled as OS X universal, but if you use NAIF's version, uncomment the following line.  i386 can be changed to x86_64 as appropriate for your system.
#ENV['ARCHFLAGS'] = '-arch i386' if RUBY_PLATFORM =~ /darwin/

require 'mkmf'

$defs.push("-std=gnu99")
$defs.push("-Wall")
$defs.push("-Werror")
$defs.push("-pedantic")

arch = RUBY_PLATFORM
arch = "darwin" if arch =~ /darwin/

gem_root = File.expand_path(File.dirname(__FILE__))

dir_config("cspice", File.join(gem_root, 'ext', 'include', 'cspice'), File.join(gem_root, 'ext', 'lib', arch))

find_header("SpiceUsr.h")
find_library("cspice", "furnsh_c")

create_makefile("spice")

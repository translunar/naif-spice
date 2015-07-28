# -*- ruby -*-

require 'rake'

desc "Load SPICE inside of a pry console"
task :pry do |task|
  cmd = [ 'pry', "-r './spice'" ]
  run *cmd
end

desc "Load SPICE inside of an irb console"
task :console do |task|
  cmd = [ 'irb', "-r './spice'" ]
  run *cmd
end

task :default => :console

def run *cmd
  sh(cmd.join(" "))
end

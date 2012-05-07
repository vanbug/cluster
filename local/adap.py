#!/usr/bin/python
import commands
fileA=open("~/chip/src/local/adapters.txt","r")

for line in fileA:
	from subprocess import call
	call(["sh","/home/genom/sukhdeeps/src/shell/yo.sh",% (line)])



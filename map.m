function out = map(in,inlow,inhigh,outlow,outhigh)
out=(in-inlow)/(inhigh-inlow)*(outhigh-outlow)+outlow;
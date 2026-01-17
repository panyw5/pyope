#!/usr/bin/env wolframscript

(* Set directory to load packages *)
SetDirectory["/Users/lelouch/pyope/.claude/skills/voa/computations"];

(* Load packages *)
Get["ope.math.m"];
Get["OPEDefs.m"];

(* Settings *)
lev = 4;

(* Output partitions *)
Print["lev = ", lev];
Print["2*lev = ", 2*lev];
part = IntegerPartitions[2*lev];
Print["IntegerPartitions[2*lev] = ", part];
Print["1/2 * IntegerPartitions[2*lev] = ", 1/2 part];
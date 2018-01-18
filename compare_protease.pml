color cyan
select protease, chain A
select peptide, chain B
color yellow, peptide
select cat_res, res 72+154
color magenta, cat_res
load a_to_s_ly104_WT.pdb
color green, a_to_s_ly104_WT
util.cnc
select des_res, res 55-58+70+73+150-153+155+170-175
select shell_res, 52+59+71+74+147-149+169+176-179+181-183
select a55, res 55
select q56, res 56
select t57, res 57
select f58, res 58
select v70, res 70
select g73, res 73
select l150, res 150
select k151, res 151
select g152, res 152
select s153, res 153
select g155, res 155
select r170, res 170
select a171, res 171
select a172, res 172
select v173, res 173
select c174, res 174
select t175, res 175
hide
show cartoon
show sticks, cat_res
show lines, peptide des_res
set line_width, 3
hide everything, elem H
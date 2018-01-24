load a_to_s_ly104_WT.pdb
color cyan
select protease, chain A
select peptide, chain B
select cat_res, res 72+96+154
select des_res, res 58+138+147+150-153+169+171-174+183
select shell_res, res 51-57+59+69-71+73-75+121+123-124+134-135+137+139+146+148-149+155-156+158+168+170+175-182
color yellow, peptide
color magenta, cat_res
color purpleblue, des_res
color green, a_to_s_ly104_WT
util.cnc

#select a55, res 55
#select q56, res 56
#select t57, res 57
#select f58, res 58
#select v70, res 70
#select g73, res 73
#select l150, res 150
#select k151, res 151
#select g152, res 152
#select s153, res 153
#select g155, res 155
#select r170, res 170
#select a171, res 171
#select a172, res 172
#select v173, res 173
#select c174, res 174
#select t175, res 175

hide
show cartoon
show sticks, cat_res 
show lines, peptide des_res
set line_width, 3
hide everything, elem H
set seq_view, 1
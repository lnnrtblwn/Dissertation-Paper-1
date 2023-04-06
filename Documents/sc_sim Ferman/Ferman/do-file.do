
clear

mat R = J(20,4,.)


import delimited "/Users/bferman/Dropbox/Projects/SC with large J and T/Simulations_Final/Table_1/Output/results.csv"

local type = "sc" // Columns 1 to 4: "sc"; columns 5 to 8: "ols_unconstrained"; columns 9 to 12: "ols_add_up"

gen mu_hat_1 = `type'_sum_w_lambda1
gen mu_hat_2 = `type'_sum_w_lambda2
gen effect = `type'_eff

preserve


keep if t0==j+5

local row=1
local col = 1
foreach x in 4 10 50 100  {

summ mu_hat_1 if j==`x'
mat R[`row',`col'] = r(mean)
mat R[`row'+1,`col'] = r(sd)

summ mu_hat_2 if j==`x'
mat R[`row'+3,`col'] = r(mean)
mat R[`row'+4,`col'] = r(sd)

summ effect if j==`x'
mat R[`row'+6,`col'] = r(sd)



local ++col 

}
restore


preserve 

keep if t0==2*j


local row=12
local col = 1
foreach x in 4 10 50 100  {

summ mu_hat_1 if j==`x'
mat R[`row',`col'] = r(mean)
mat R[`row'+1,`col'] = r(sd)

summ mu_hat_2 if j==`x'
mat R[`row'+3,`col'] = r(mean)
mat R[`row'+4,`col'] = r(sd)

summ effect if j==`x'
mat R[`row'+6,`col'] = r(sd)



local ++col 

}

restore


mat li R




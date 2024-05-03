
some_function(2, option = "r")

some_function(seq(200), option = "eigen")

# checking epi spread function from finalsize
polymod <- socialmixr::polymod
contact_matrix <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40)
)
demography <- contact_matrix$participants$proportion
demography = demography / sum(demography)
contact_matrix <- t(contact_matrix$matrix)
# scale contact matrix cols by prop of dem
for(i in seq(ncol(contact_matrix))) contact_matrix[, i] = contact_matrix[, i] / demography[i]

contact_matrix = contact_matrix / (max(Re(eigen(contact_matrix)$values)))

p_susceptibility <- matrix(1, 3, 3)
susceptibility = matrix(1, 3, 3)

devtools::load_all()

solve_final_size_internal(
    contact_matrix = contact_matrix,
    demography = demography,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
)

# Function - get species sci name if not provided, includes nei species
get_names <- function(x) {
  subfam <- species_list(Subfamily = x)
  fam <- species_list(Family = x)
  gen <- species_list(Genus = x)
  cls <- species_list(Class = x)
  ord <- species_list(Order = x)

  final <- c(subfam, fam, gen, cls, ord)

  if(length(final) > 0) {
  return(final)
  }
}

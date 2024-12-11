#####################################################
# ~ cutter: normaliseCoverage ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

# function is to re-create "mut" table
# but as if every sample had exactly the same coverage
# it is useful for histograms, for example
# it uses the 'freq' column to calculate predicted counts
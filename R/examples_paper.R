#**********************************************************#
#********************* CREATE TABLES **********************#
#**********************************************************#

# ------------------------ Afghan Girls RCT
afghan_obj <- list(data = afghan,
                   y_name = 'testscore',
                   T_name = 'enrolled',
                   z_name = 'buildschool',
                   controls = c('headchild', 'age', 'yrsvill', 'farsi', 'tajik',
                                'farmers', 'agehead', 'educhead', 'nhh', 'land',
                                'sheep', 'distschool', 'chagcharan'),
                   evaluateInterior=FALSE,
                   example_name = 'Afgan Girls RCT',
                   option = c(NA,NA,"PequalPstar"),
                   a0bound = c(1,1,1),
                   a1bound = c(1,1,1),
                   dTstar_tilde_range = matrix(c(NA,NA,
                                                 0,1,
                                                 0,1),ncol=2,byrow=TRUE))

example_afghan <- do.call(makeExample,afghan_obj)

# ------------------------ Smoking and BMI
smoking_obj <- list(data = smoking,
                    y_name = 'BMI',
                    T_name = 'quit',
                    z_name = 'program',
                    evaluateInterior=TRUE,
                    example_name = 'Smoking \\& BMI',
                    option = c(NA,NA),
                    a0bound = c(1,1),
                    a1bound = c(0,0),
                    dTstar_tilde_range = matrix(c(NA,NA,
                                                  -1.5,0),ncol=2,byrow=TRUE))

example_smoking <- do.call(makeExample,smoking_obj)

setwd('C:/Users/mmm/Dropbox/UPenn/Research Assistant/Frank & Camilo/Bounds_binary case')

makeTable(file="tables_binary_case.tex",
          example_afghan, example_smoking)
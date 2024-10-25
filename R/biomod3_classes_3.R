
## -------------------------------------------------------------------------- #
## 1. MS.models.out -----------------------------------------------------
## -------------------------------------------------------------------------- #

setClass("MS.models.out",
         representation(ms.project = 'character',
                        modeling.id = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        data.type = 'character',
                        expl.var.names = 'character',
                        summary.models.computed = 'ANY',
                        has.evaluation.data = 'logical',
                        scale.models = 'logical'),
         prototype(modeling.id = as.character(format(Sys.time(), "%s")),
                   dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   summary.models.computed = NULL,
                   has.evaluation.data = FALSE,
                   scale.models = TRUE),
         validity = function(object){ return(TRUE) } )

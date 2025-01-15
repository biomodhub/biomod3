
## biomod3_internal ------------
utils::globalVariables(names = c("categories", "cats", "sub.i", "melt", "i", "i.dim1", "i.dim2", "i.dim3"
                                 ))

## BIOMOD_ProjectionWrap -------
utils::globalVariables(names = c("head", "new", "models.chosen.single", "models.chosen.ens"))
                       
## BIOMOD_Wrap -----------------
utils::globalVariables(names = c("bm.format", "expl.var"))
                       
## MS_EnsembleForcasting -------
utils::globalVariables(names = c("sp"))
                       
## MS_EnsembleModeling ---------
## MS_FormatingData ------------
## MS_Modeling -----------------
## MS_Projection ---------------
## bm_SpeciesParameters --------
utils::globalVariables(names = c("user.table"))

## BIOMOD.wrap.out -------------
## MS.ensemble.models.out ------
## MS.formated.data ------------
utils::globalVariables(names = c("has.mask", "has.mask.eval", "this_dataset", "value",
                                 "y", "resp"))

## MS.models.out ---------------
utils::globalVariables(names = c("mod"))

## MS.projection.out -----------
utils::globalVariables(names = c("proj", "pred"))
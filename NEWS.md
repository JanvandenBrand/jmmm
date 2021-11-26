# jmmm 0.6

## Major changes
Exported the following functions:
* `get_mmm_start_values()`
* `get_mmm_derivatives()`

Added `average_pairwise_models()`

Updated `get_start_values()` to no longer stop execution of `mmm_model()` if fitting a model fails. Instead a placeholder message is included.

Together these changes should help debug and update the multivariate mixed model if errors occur during model fitting. For an example on how to use these functions for debugging see the README file.

# jmmm 0.5

## Minor changes
Export family function

# jmmm 0.4

## Minor changes
Added error handling to retrieve the model on which fitting failed when using a parallel back end.

# jmmm 0.3

## Minor changes
Updated handling of `progressr` to show more accurate status updates on fitting pairwise models in parallel. 

# jmmm 0.2

## Major changes
Update to the mmm_model() function. It no longer requires the user to specify the parallel back end as an argument. Instead the user should specificy their preferred back end prior to running the function. This way users have more control over parallelization in their specific setting.

## Minor changes
* Added examples to documentation.
* Added versions to package dependencies. 

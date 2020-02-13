setClass(
    'GenomicSignature',
    representation = representation(
        meta = 'list',
        components = 'matrix',
        "VIRTUAL"
    ),
    prototype = list(
        meta = list(),
        components = matrix()
    )
)

setClass(
    'PCAGenomicSignature',
    contains = c('GenomicSignature')
)

setClass(
    'ClusteredPCAGenomicSignature',
    contains = 'GenomicSignature',
    slots = c(
        k = 'integer'
    ),
    prototype = list(
        k = 0L,
        meta = list(),
        components=matrix()
    )
)
    

###################################
## Working with component matrix ##
###################################

# accessor

setGeneric(
    'components',
    def = function(object) {
        standardGeneric('components')
    }
)

setMethod(
    'components',
    signature = "GenomicSignature",
    definition = function(object) {
        return(object@components)
    }
)

# setter

setGeneric(
    'components<-',
    def = function(object, value) {
        standardGeneric('components<-')
    }
)


setMethod(
    'components<-',
    signature = c("GenomicSignature"),
    function(object, value) {
        object@components <- value
        return(object)
    }
)

#############################
## Getting "k" from object ##
#############################

setGeneric(
    'clusterCount',
    def = function(object) {
        standardGeneric('clusterCount')
    }
)

setMethod(
    'clusterCount',
    signature = "ClusteredPCAGenomicSignature",
    definition = function(object) {
        return(object@k)
    }
)

## z = new('PCAGenomicSignature')
## components(z)
## components(z) = matrix(0, nrow=10, ncol=10)
## z2 = new('ClusteredPCAGenomicSignature')
## components(z2)
## clusterCount(z2)



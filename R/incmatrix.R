#' @import dplyr
#' @import sf
#' @export
#' @title
#' Function that build incidence matrix!
#' @description
#' The binary relation that represents the incidence matrix is 'Polygon A shares a border with Polygon B.
#' It returns a matrix with 0 and 1 that specify if there is incidency between poligons
#' @details
#' ESPECIFICAR TODO LOS DETALLES TËCNICOS NECESARIOS
#' PUEDEN USARSE VARIAS LINEAS E INTRODUCIR ECUACIONES
#' \deqn{\hat{\sigma}^2_y = \frac{\sum_s w_k (y_k - \bar{y}_w)^2}{\sum_s w_k}}
#' @author Javier Nuñez <Javier_Nunez at inec.gob.ec>
#' @author Christian Garces <Christian_Garces at inec.gob.ec>
#' @param pol Non-adjacent polygons generally refer to blocks.
#' @param tol Tolerance in meters.
#' @param id ID of the original Primary Sampling Units
#'
#' @references
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas.
#' Valliant, R, et. al. (2013), \emph{Practical tools for Design and Weighting Survey Samples}. Springer
#' @return
#' @export
#'
#' @examples

incmatrix <- function(pol, tol=10, id = NULL){
  aux <- pol %>%
    rename(id = {{id}})
  aux1 <- st_intersection(aux,aux) %>%
    group_by(id) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    filter(!(id==id.1 & n > 1)) %>%
    mutate(largo=as.numeric(st_length(.)),
           frecuencia = 1) %>%
    filter(largo > tol | largo == 0) %>%
    as.data.frame() %>%
    select(id, id.1, frecuencia)

  o <- aux1 %>%
    arrange(id.1) %>%
    pivot_wider(names_from = id.1, values_from = frecuencia) %>%
    arrange(id) %>%
    select(id, .$id) %>%
    as.data.frame()

  q <- data.matrix(select(o, -id)) %>%
    replace(is.na(.), 0)

  rownames(q) <- colnames(q)

  return(q)
}

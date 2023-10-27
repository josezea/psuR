#' @import dplyr
#' @import sf
#' @export
#' @title
#' Function that extends the polygons based on cartographic information
#' @description
#' This function allows for extending the polygons prior to the formation of Primary Sampling Units
#' @details
#' ESPECIFICAR TODO LOS DETALLES TËCNICOS NECESARIOS
#' PUEDEN USARSE VARIAS LINEAS E INTRODUCIR ECUACIONES
#' \deqn{\hat{\sigma}^2_y = \frac{\sum_s w_k (y_k - \bar{y}_w)^2}{\sum_s w_k}}
#' @author Javier Nuñez <Javier_Nunez at inec.gob.ec>
#' @author Christian Garces <Christian_Garces at inec.gob.ec>
#' @param poli Non-adjacent polygons generally refer to blocks.
#' @param profile The polygon that defines the area where the extension work will take place
#' @param id ID of the original Primary Sampling Units
#' @param boundary It's a multi-polygon that represents rivers or boundaries, used to avoid connections between blocks separated by them, and it can or cannot not be contained within the profile
#' @param buf Buffer Width is a scalar value that determine the quality of the polygon extension. Higher quality requires more processing
#' @param density TDensity buffer is a scalar value that determine the quality of the polygon extension. Higher quality requires more processing
#'
#' @references
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas.
#' Valliant, R, et. al. (2013), \emph{Practical tools for Design and Weighting Survey Samples}. Springer
#' @return
#' @export
#'
#' @examples

extpol <- function(poli, profile, id = NULL, boundary = NULL, buf = 5, density = 0.1){
  # source("rutinas/funciones/box.R")

  box <- function(puntos){
    bb <- st_bbox(puntos) #  sf::st_bbox
    p<-matrix(
      c(bb["xmin"],bb["ymin"],
        bb["xmin"],bb["ymax"],
        bb["xmax"],bb["ymax"],
        bb["xmax"],bb["ymin"],
        bb["xmin"],bb["ymin"]),
      ncol=2,byrow = T
    )
    box <- st_polygon(list(p)) # sf::st_polygon
    return(box)
  }

  print(paste0("The function has started ", Sys.time()))

  poli <- poli %>%
    rename(id = {{ id }})

  index <- unique(poli$id)
  quntosi <- vector("list", 0)

  for(i in 1:length(index)){

    pol <- poli %>%
      filter(id == index[i])

    pol1 <- st_buffer(pol, -buf) %>%
      summarise()

    pol2 <- st_difference(pol, pol1)

    quntos1 <- pol2 %>%
      st_simplify(dTolerance = 0.9) %>%
      mutate(area = as.numeric(st_area(.)),
             npuntos = ceiling(area/(1/density))) %>%
      st_cast("MULTIPOLYGON") %>%
      st_sample(size = sum(.$npuntos), type = "hexagonal") %>%
      st_as_sf() %>%
      st_cast("POINT") %>%
      st_coordinates() %>%
      as.data.frame() %>%
      st_as_sf(coords = c("X","Y"))

    st_crs(quntos1) <- st_crs(pol)

    quntos_iden <- quntos1 %>%
      st_join(pol %>% select(id), join = st_within) %>%
      filter(!is.na(id))

    names(quntos_iden)[names(quntos_iden) == attr(quntos_iden, "sf_column")] = attr(pol, "sf_column")
    st_geometry(quntos_iden) <- attr(pol, "sf_column")

    quntosi[[i]] <- pol %>%
      #st_buffer(-0.2) %>% #no estaba comentada
      st_simplify(dTolerance = 0.9) %>%
      st_as_sf() %>%
      st_cast("MULTIPOINT") %>%
      st_cast("POINT") %>%
      filter(!is.na(id)) %>%
      select(id) %>%
      rbind(quntos_iden)

    print(length(index) - i)

  }

  quntos <- do.call(rbind, quntosi)

  puntos <- st_union(quntos)

  st_crs(puntos) <- st_crs(poli)
  st_crs(quntos) <- st_crs(poli)
  print(paste0("The points have been created: ", Sys.time()))
  #Calculamos los polígonos de voronoi
  caj <- st_sfc(box(puntos), crs = st_crs(poli))
  voronoi <- st_voronoi(puntos, caj) %>%
    st_cast() %>%
    st_as_sf()
  voronoi <- voronoi %>%
    # data.frame(geometry = .) %>%
    # st_sf(.) %>%
    st_join(.,quntos, join = st_contains) %>%
    #calculo de área
    mutate(area = as.numeric(st_area(.)))
  #Corrección voronoi
  # for (i in 1:(dim(voronoi)[1])){
  #   if(is.na(voronoi$man[i])){
  #     voronoi$man[i] <- voronoi$man[i-1]
  #   }
  # }
  print(paste0("The Voronoi polygons have been created: ", Sys.time()))
  vorpun <- voronoi
  # disolver por manzana
  disolver <- voronoi %>%
    st_buffer(0) %>% #1
    st_make_valid() %>%
    group_by(id) %>%
    summarise(np=n()) %>%
    st_buffer(0) %>% #-1
    st_cast() %>%
    filter(!is.na(id))
  st_crs(disolver) <- st_crs(pol)
  print(paste0("The Voronoi polygons have been dissolved: ", Sys.time()))
  # se intersecta con el profile
  polext <- st_intersection(disolver,profile) %>%
    group_by(id) %>%
    summarise() %>%
    st_cast()
  print(paste0("The polygons have been cut by the profile: ", Sys.time()))
  # polext <- polext %>%
  #   mutate(loli = st_is_valid(.))


  if(!is.null(boundary)){
    polext <- polext %>%
      st_difference(boundary) %>%
      st_buffer(0) %>%
      #Validación de geometría
      st_make_valid() %>%
      #Trasnformación a multipoligono
      st_cast("MULTIPOLYGON") %>%
      #Transformación a poligono
      st_cast("POLYGON") %>%
      #Emparejamiento geográfico para identicar pedazos de manzanas a utilizar
      st_join(poli) %>%
      filter(id.x == id.y) %>%
      #Disolver a nivel de manzana
      group_by(id = id.x) %>%
      summarise() %>%
      st_make_valid() %>%
      st_cast("MULTIPOLYGON")
  }

  polext <- polext %>%
    rename({{ id }} := id) %>%
    st_as_sf()

  names(polext)[names(polext)==attr(polext, "sf_column")] = "geom"
  st_geometry(polext)="geom"
  print(paste0("The intermediate dataframes have been generated: ", Sys.time()))

  List <- list(polext, quntos)

  rm(quntosi, pol, pol1, pol2, quntos1, quntos_iden, puntos, voronoi, vorpun, disolver, quntos, polext)

  return(List)
}

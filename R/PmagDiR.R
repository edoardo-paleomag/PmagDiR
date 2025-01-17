#temporary WebDiR NOT READY, DO NOT USE
WebDiR <- function(){
  library(PmagDiR)
  library(plyr)
  library(dplyr)
  library(shiny)
  library(shinyWidgets)
  library(DT)
  library(shinyhelper)


  ui <- fluidPage(
    #navbarPage("Web-DiR, by Edoardo Dallanave"),
    tabsetPanel(
      type = "tabs",
      tabPanel("Directions analysis",
               tabsetPanel(
                 tabPanel("Vector end-points interpolation",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                                         fluidRow(
                                           column(6,fileInput("All_Zijd",label = "Load demag data",multiple = T)),
                                           column(6,selectInput("Zijd_f_type",label = "File type",
                                                                choices = list("WebDiR"=1,"Lamont"=2,"Bremen(.cor)"=3,"IODP_JR6A_Expanded"=4,"CIT_samples"=5),selected = 1))
                                         ),
                                         fluidRow(
                                           column(6,selectInput("Zijd_Stereo_shift",label = "Diagram",
                                                                choices = list("N-Right"=1,"N-Up"=2, "Equal area"=3), selected = 1)),
                                           column(6,selectInput(inputId = "VEPcoordinates",label = "Coordinates",
                                                                choices = list("Specimen"=1,"Geographic"=2,"Tilt Corr."=3),selected = 3))
                                         ),
                                         fluidRow(
                                           column(6,textInput(inputId = "textunit",label = "Unit",value = "A/m")),
                                           column(6,selectInput("VEPticks",label = "Ticks",
                                                                choices = list("x0.05"=1,"x0.1"=2,"x0.25"=3,"x0.5"=4,"x1.0"=5,"No ticks"=6),selected = 4))
                                         ),
                                         fluidRow(
                                           column(6,selectInput("anchor",label = "Interpolation",
                                                                choices = list("PCA Free"=1,"PCA Anch."=2,"PCA Or. Incl."=3, "Fisher"=4, "G. Circle"=5),selected = 1)),
                                           column(6, textInput("comp_name",label = "Component name",value = "Ch"))
                                         ),
                                         fluidRow(
                                           column(12,h4("List of loaded specimens"))
                                         ),
                                         fluidRow(
                                           column(12,DT::dataTableOutput("samples_list"))
                                         )
                            ),
                            mainPanel(
                              fluidRow(
                                actionButton(inputId = "del_VEPs",label = "Delete selection"),
                                actionButton(inputId = "restore_VEPs",label = "Restore data"),
                                actionButton(inputId = "runVEPstat",label = "Run stat"),
                                actionButton(inputId = "save_PCA",label = "Save stat"),
                                downloadButton("export_PCA",label = "Export all saved directions"),
                                downloadButton("export_Zijd",label = "Export figure")
                              ),
                              fluidRow(column(3,textOutput("Zijd_Unit")),
                                       column(9,textOutput("PCA_result"))
                              ),
                              column(11,plotOutput("zijderveld",brush = brushOpts(id = "plot_brush", fill = NA))),
                              column(1, DT::dataTableOutput("sampledat",width = 100))
                            )
                          )
                 ),
                 tabPanel("All saved samples",
                          sidebarLayout(
                            sidebarPanel(width = 6,
                                         fluidRow(
                                           column(4,fileInput(inputId = "import_PCA",label = "Import saved directions")%>%
                                                    helper(type = "inline",
                                                           title = "Format file",
                                                           content = c(
                                                             "File as exported from Vector end-points interpolation page"),
                                                           size = "m",fade = T)),
                                           column(4,textInput(inputId = "sel_interpol_name",label = "Name of exported file",value = "Directions")),
                                           column(4,selectInput(inputId = "EAcoordinates",label = "Coordinates",
                                                                choices = list("Geographic"=1,"Tilt Corr."=2),selected = 2))
                                         ),
                                         fluidRow(DT::dataTableOutput(outputId = "saved_interpol"))
                            ),
                            mainPanel(width = 6,
                                      fluidRow(actionButton(inputId = "del_interpol",label = "Delete selection"),
                                               actionButton(inputId = "undel_interpol",label = "undo delete"),
                                               downloadButton("export_interpol",label = "Export sel. directions"),
                                               actionButton(inputId = "comb_DI_GC",label = "Combine DI & GC"),
                                               actionButton(inputId = "save_GC",label = "Add GC dirs to list"),
                                               actionButton(inputId = "GC_erase",label = "Clear GC dirs from plot")),
                                      plotOutput("saved_interpol_EA"))
                          )
                 )
               )
      ),
      tabPanel("Directions display and filter",
               sidebarLayout(
                 sidebarPanel(
                   fluidRow(
                     column(6,fileInput("file", label= "Directions file input")),
                     column(6,textInput("fileN",label = "Export name",value = "Site"))
                   ),
                   fluidRow(
                     column(8,selectInput("filetype", label = "DIRECTIONS FILE TYPE",
                                          choices = list("Dec, Inc "=1,"G_dec, G_inc, B_az, B_plung"=2,"G_dec, G_inc, B_dec, B_inc"=3,
                                                         "Web-DiR"=4),selected = 1) %>%
                              helper(type = "inline",
                                     title = "Format file",
                                     content = c(
                                       "HELP TO BE FIXED",
                                       "Dec, Inc (, Bed_az, Bed_plung)= declination and inclination, with optional bedding azimuth and plunge
                                               if data are in geographic coordinates or for the 'unstraining' process.",
                                       "",
                                       "G_dec, G_inc, B_dec, B_inc= Geographic coordinates declination and inclination,
                                               followed by bedding corrected declination and inclination.",
                                       "",
                                       "An optional third (or fifth) column containing the stratigraphic position in meters will be used for
                                              the magnetic polarity stratigraphy plot. Any other column will be ignored."),
                                     size = "m",fade = T)
                     ),
                     column(4,selectInput("coord", label= "Coordinates",
                                          choices = list("Geographic"=1, "Bedding"=2,"Specimen"=3),selected = 1))
                   ),
                   fluidRow(
                     column(6,numericInput("lat",label="Site latitude",value=0)),
                     column(6,numericInput("long",label="Site longitude",value=0))
                   ),
                   fluidRow(
                     column(6,selectInput("fisher", label= "Mean",
                                          choices = list("None"=1, "Fisher" = 2, "Elliptic" = 3,
                                                         "Inc. only single mode"=4,"Inc. only bimodal"=5,
                                                         "Arithm. single mode"=6, "Arithm. bimodal"=7),selected = 1)),
                     column(6,selectInput("mode", label= "Mode",
                                          choices = list("Bimodal"=1, "Mode 1" = 2, "Mode 2" = 3, "All down"=4, "All up"=5),selected = 1))
                   ),
                   fluidRow(
                     column(4,selectInput("colD", label= "Color down",
                                          choices = list("black"=1, "blue" = 2, "red" = 3,"green"=4),selected = 2)),
                     column(4,selectInput("colU", label= "Color up",
                                          choices = list("white"=1, "cyan" = 2, "pink" = 3,"light green"=4),selected = 2)),
                     column(4,selectInput("sym", label= "Symbol",
                                          choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1))
                   ),
                   fluidRow(
                     column(6,selectInput("cutoff", label= "Cut-off",
                                          choices = list("None"=1, "Vandamme dynamic" = 2, "Vandamme static" = 3,
                                                         "Fixed dynamic"=4, "Fixed static"=5,
                                                         "Cut up-pointing"=6,"Cut down-pointing"=7,"Cut fixed inc."=8),selected = 1)),
                     column(6, numericInput("VGP_fixed", label= "VGP fixed-filter radius", value=45)),
                   ),
                   fluidRow(
                     column(6, numericInput("MinInc",label = "Min Inc. filt.",value = 0)),
                     column(6, numericInput("MaxInc",label = "Max Inc. filt.",value = 0)),
                   ),
                   fluidRow(
                     tableOutput("stats")
                   ),
                   fluidRow(
                     h4(textOutput("inc_warn"))
                   ),
                   br(),
                   fluidRow(
                     column(12,actionButton("resetDir",label = "Delete input file", width = "100%"))
                   ),
                   br(),
                 ),
                 mainPanel(
                   fluidRow(downloadButton("exportG","Export graph"),
                            downloadButton("exportS","Export stat"),
                            downloadButton("exportDI","Export directions")),
                   column(1),
                   plotOutput("directions")
                 )
               )
      ),
      tabPanel("Reversal test",
               sidebarLayout(
                 sidebarPanel(width = 3,
                              fluidRow(
                                column(6,textInput("fileN_RT",label = "Export name",value = "Site")),
                                column(6,numericInput("revnb", label="Bootstraps n.",value=1000))),
                              fluidRow(
                                column(12,actionButton("revgo", label= "Execute",width = "100%"))),
                              fluidRow(
                                br(),
                                column(12,progressBar(
                                  id = "mode1",
                                  value = 0,total=1000,
                                  title = "Bootstrap Mode 1",
                                  display_pct = TRUE))),
                              fluidRow(
                                br(),
                                column(12,progressBar(
                                  id = "mode2",
                                  value = 0, total=1000,
                                  title = "Bootstrap Mode 2",
                                  display_pct = TRUE
                                ))
                              ),fluidRow(
                                tableOutput("revstat")
                              )
                 ),
                 mainPanel(
                   fluidRow(downloadButton("revexpG","Export graph"),
                            downloadButton("revexpS","Export stat")),
                   plotOutput("revtest")
                 )
               )
      ),
      tabPanel("Distribution shape",
               sidebarLayout(
                 sidebarPanel(width = 3,
                              fluidRow(
                                column(6,textInput("fileN_EI",label = "Export name",value = "Site")),
                                column(6,numericInput("EIbootnb", label="Bootstraps n.",value=5000))),
                              fluidRow(
                                column(6,selectInput("EIyesnoboot", label="Bootstrap?",
                                                     choices = list("No"= 1, "Yes" = 2),selected = 2)),
                                column(6,numericInput("EIconf", label="confidence %",value = 95))),
                              fluidRow(
                                column(12,actionButton("EIbootgo", label= "Execute",width = "100%"))),
                              fluidRow(
                                br(),
                                column(12,progressBar(
                                  id = "EIbootstrap",
                                  value = 0,total=1000,
                                  title = "Bootstrap",
                                  display_pct = TRUE
                                ))
                              ),
                              br(),
                              fluidRow(
                                column(6, downloadButton("EIbootG","Export graph")),
                                column(6, downloadButton("EIbootS","Export stat"))),
                              br(),
                              br(),
                              fluidRow(
                                tableOutput("EIbootStat")
                              )
                 ),
                 mainPanel(
                   column(1),
                   plotOutput("EIboot")
                 )
               )
      ),
      tabPanel("Inclination flattening",
               tabsetPanel(
                 tabPanel("EI graph",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                                         fluidRow(
                                           column(6,textInput("fileN_FF",label = "Export name",value = "Site")),
                                           column(6,numericInput("ffindboot", label="Bootstraps number",value=1000))
                                         ),
                                         fluidRow(
                                           column(6,selectInput("ffindyesnoboot", label="Bootstrap?",
                                                                choices = list("No"= 1, "Yes" = 2),selected = 1))
                                         ),
                                         fluidRow(
                                           column(12,actionButton("ffindgo", label= "Execute",width = "100%"))
                                         ),
                                         fluidRow(
                                           br(),
                                           column(12,progressBar(
                                             id = "ffindbootstrap",
                                             value = 0,total=1000,
                                             title = "Valid bootstraps",
                                             display_pct = TRUE
                                           ))
                                         ),
                                         fluidRow(h5(textOutput("validboots"))),
                                         br(),
                                         fluidRow(
                                           column(6, downloadButton("ffindG","Export graph")),
                                           column(6, downloadButton("ffindS","Export stat"))
                                         ),
                                         br(),
                                         fluidRow(
                                           tableOutput("ffindStat")
                                         )
                            ),
                            mainPanel(
                              column(1),
                              plotOutput("ffindgraph")
                            )
                          )
                 ),
                 tabPanel("Unflattened directions",
                          sidebarLayout(
                            sidebarPanel(
                              fluidRow(
                                column(8,selectInput("fisher_FF", label="Mean",
                                                     choices = list("None"=1, "Fisher" = 2, "Elliptic" = 3),selected = 1)),
                                column(4,selectInput("mode_FF", label= "Mode",
                                                     choices = list("Bimodal"=1, "Mode 1" = 2, "Mode 2" = 3, "All down"=4, "All up"=5),selected = 1))
                              ),
                              fluidRow(
                                column(4,selectInput("colD_FF", label= "Color down",
                                                     choices = list("black"=1, "blue" = 2, "red" = 3,"green"=4),selected = 2)),
                                column(4,selectInput("colU_FF", label= "Color up",
                                                     choices = list("white"=1, "cyan" = 2, "pink" = 3,"light green"=4),selected = 2)),
                                column(4,selectInput("sym_FF", label= "Symbol",
                                                     choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1))
                              ),
                              tableOutput("stats_FF")
                            ),
                            mainPanel(
                              fluidRow(downloadButton("ffindG_FF","Export graph"),
                                       downloadButton("ffindS_FF","Export stat"),
                                       downloadButton("ffindDI_FF","Export directions")),
                              column(1),
                              plotOutput("directions_FF")
                            )
                          )
                 )
               )
      ),
      tabPanel("Strain removal",
               tabsetPanel(
                 tabPanel("Initial evaluation",
                          sidebarLayout(
                            sidebarPanel(
                              fluidRow(
                                column(6,textInput("fileN_str",label = "Export name",value = "Site")),
                                column(6,fileInput("str_mtrx",label = "Strain matrix input"))
                              ),
                              tableOutput("straindirs"),
                              fluidRow(
                                column(4,selectInput("str_m_type",label = "Type",
                                                     choices = list("Vi,dec,inc"=1,"Kij"=2, "3x3"=3),selected = 1)),
                                column(4,selectInput("brk",label = "Break",
                                                     choices = list("None"=1,"Cross TK03"=2,"Max Edec"=3,"Min Edec"=4),selected = 1)),
                                column(4,numericInput("increm",label = "Increments n.",value = 50))
                              ),
                              fluidRow(
                                column(4,numericInput("lin",label = "Target lineation",value = 1)),
                                column(4,numericInput("fol",label = "Target foliation", value = 1))
                              ),
                              br(),
                              fluidRow(
                                column(12,progressBar(
                                  id = "nincrements",
                                  value = 0,total=50,
                                  title = "Total 'unstrain' increments",
                                  display_pct = TRUE
                                ))
                              ),
                              textOutput("finalpars"),
                              br(),
                              fluidRow(
                                column(6,actionButton("unstr_GO",label = "Execute",width = "100%")),
                                column(6,downloadButton("unstrG",label = "Export graph"))
                              )
                            ),
                            mainPanel(
                              column(1),
                              plotOutput("unstrainplot")
                            )
                          )
                 ),
                 tabPanel("Bootstrap analysis",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                                         fluidRow(
                                           column(6,textInput("fileN_strBoot",label = "Export name",value = "Site")),
                                           column(6,numericInput("strboot", label="Boot. number",value=1000))
                                         ),
                                         fluidRow(
                                           column(6,numericInput("incremboot",label = "Increments n.",value = 20)),
                                           column(6,selectInput("strainhist",label = "Histogram?",
                                                                choices = list("Yes"=1,"No"=2),selected = 1))
                                         ),
                                         br(),
                                         fluidRow(
                                           column(12,progressBar(
                                             id = "str_boots",
                                             value = 0,total=1000,
                                             title = "Bootstrapped pseudosamples",
                                             display_pct = TRUE)
                                           )
                                         ),
                                         h5("Bootstrap process is time demanding"),
                                         br(),
                                         fluidRow(
                                           column(6,actionButton("unstr_boot_GO",label = "Execute",width = "100%"))
                                         ),
                                         br(),
                                         fluidRow(
                                           column(6, downloadButton("unstrbootG", label = "Export Graph")),
                                           column(6, downloadButton("unstr_boot_result", label = "Export Stat"))
                                         ),
                                         br(),
                                         fluidRow(
                                           tableOutput("unstr_boot_stat")
                                         )
                            ),
                            mainPanel(
                              column(1),
                              plotOutput("unstrainboot")
                            )
                          )
                 ),
                 tabPanel("Corrected directions",
                          sidebarLayout(
                            sidebarPanel(
                              fluidRow(
                                column(4,selectInput("fisher_str", label="Mean",
                                                     choices = list("None"=1, "Fisher" = 2, "Elliptic" = 3),selected = 1)),
                                column(4,selectInput("mode_str", label= "Mode",
                                                     choices = list("Bimodal"=1, "Mode 1" = 2, "Mode 2" = 3, "All down"=4, "All up"=5),selected = 1)),
                                column(4,selectInput("coord_str", label = "Coordinates",
                                                     choices = list("Geographic"=1,"Bedding"=2), selected = 1))
                              ),
                              fluidRow(
                                column(6,selectInput("cutoff_str", label= "Cut-off",
                                                     choices = list("None"=1, "Vandamme dynamic" = 2, "Vandamme static" = 3,
                                                                    "Fixed dynamic"=4, "Fixed static"=5),selected = 1)),
                                column(6, numericInput("VGP_fixed_str", label= "VGP fixed-filter radius", value=45)),
                              ),
                              fluidRow(
                                column(4,selectInput("colD_str", label= "Color down",
                                                     choices = list("black"=1, "blue" = 2, "red" = 3,"green"=4),selected = 2)),
                                column(4,selectInput("colU_str", label= "Color up",
                                                     choices = list("white"=1, "cyan" = 2, "pink" = 3,"light green"=4),selected = 2)),
                                column(4,selectInput("sym_str", label= "Symbol",
                                                     choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1))
                              ),
                              tableOutput("unstr_DI_stat")
                            ),
                            mainPanel(
                              fluidRow(downloadButton("unstrdirsG","Export graph"),
                                       downloadButton("unstrdirsS","Export stat"),
                                       downloadButton("unstrdirsD","Export directions")),
                              column(1),
                              plotOutput("unstrDirs")
                            )
                          )
                 )
               )
      ),
      tabPanel("VGPs",
               tabsetPanel(
                 tabPanel("VGPs plot and average",
                          sidebarLayout(
                            sidebarPanel(width = 5,
                                         fluidRow(
                                           column(12,h4("Loaded directions  VGPs")),
                                         ),
                                         fluidRow(
                                           column(6,textInput("fileN_VGP",label = "Export name",value = "Site")),
                                           column(6,selectInput("dirs_vgp", label="Directions",
                                                                choices = list("Original"= 1,"Unflattened"= 2), selected = 1))
                                         ),
                                         fluidRow(
                                           column(4,selectInput("centercoord",label="Center of plot",
                                                                choices=list("Automatic"= 1,"Manual"= 2), selected=1)),
                                           column(4,numericInput("centerlat",label="Center lat", value=90)),
                                           column(4,numericInput("centerlong",label="Center Long", value=0))
                                         ),
                                         fluidRow(
                                           column(4,selectInput("coastyesno", label = "Plot coastline",
                                                                choices = list("Yes"=1,"no"=2),selected=1)),
                                           column(4,selectInput("vgpscolor", label= "VGP color",
                                                                choices= list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,"red"=7,"yellow"=8,"cyan"=9, "gray"=10,"white"= 11), selected=3)),
                                           column(4,selectInput("vgpssymbol",label = "VGP symbol",
                                                                choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1))
                                         ),
                                         fluidRow(
                                           column(4,selectInput("plotA95", label="Statistic",
                                                                choices=list("None"=1,"Fisher"=2,"Bootstrap"=3),selected=1)),
                                           column(4,numericInput("vgpbootn",label="Bootstraps n.",value=2000)),
                                           column(4,selectInput("VGPtype", label = "VGPs to export",
                                                                choices = list("Single mode"=1,"Bimodal"=2,"Rotated"=3),selected = 1))
                                         ),
                                         fluidRow(
                                           br(),
                                           column(12,progressBar(
                                             id = "vgpboot",
                                             value = 0,total=1000,
                                             title = "VGPs bootstrap",
                                             display_pct = TRUE))
                                         ),
                                         tableOutput("fishpole"),
                                         br(),
                                         h4(textOutput("geowarning")),
                                         h4(textOutput("coordwarning")),
                                         fluidRow(
                                           column(12,actionButton("saveVGP",label = "Add to global list",width = "100%"))
                                         ),
                                         br(),
                                         h4("Preliminary list of loaded VGPs:"),
                                         h5("(Interactive list on Multiple VGPs analysis)"),
                                         br(),
                                         tableOutput("Ext_VGP_list3")
                            ),
                            mainPanel(width = 7,
                                      fluidRow(
                                        downloadButton("VGPs_S","Export Pole"),
                                        downloadButton("VGPs_Exp","Export VGPs"),
                                        downloadButton("VGPs_G","Export graph")
                                      ),
                                      column(1),
                                      plotOutput("VGPplot")
                            )
                          )
                 ),
                 tabPanel("External VGP",
                          sidebarLayout(
                            sidebarPanel(width = 5,
                                         fluidRow(
                                           column(12,h4("Externally loaded VGPs")),
                                         ),
                                         fluidRow(
                                           column(6,fileInput("vgpfile", label = "Load VGPs file")),
                                           column(6,textInput("VGP_ext_sitename", label = "VGPs name",value = "Locality"))
                                         ),
                                         fluidRow(
                                           column(3,selectInput("VGP_ext_center",label="Center of plot",
                                                                choices=list("Automatic"= 1,"Manual"= 2), selected=1)),
                                           column(3,numericInput("VGP_ext_clat",label="Center lat", value=90)),
                                           column(3,numericInput("VGP_ext_clong",label="Center Long", value=0)),
                                           column(3,selectInput("VGP_ext_coast", label = "Coastline",
                                                                choices = list("Yes"=1,"no"=2),selected=1)),
                                         ),
                                         fluidRow(
                                           column(3,selectInput("MVGPcolor", label= "VGP color",
                                                                choices= list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,"red"=7,"yellow"=8,"cyan"=9,"gray"=10, "white"=11), selected=3)),
                                           column(3,selectInput("MVGPsymbol", label= "VGP symbol",
                                                                choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1)),
                                           column(3,selectInput("MVGP_stat", label="Statistic",
                                                                choices=list("None"=1,"Fisher"=2,"Bootstrap"=3),selected=1)),
                                           column(3,numericInput("MVGPnb_ext", label = "Bootstrap n.", value = 2000))
                                         ),
                                         fluidRow(
                                           column(4,selectInput("VGP_ext_cut_type",label = "Cut-off",choices = list("None"=1,"Vandamme"=2,"Fixed"=3),selected = 1)),
                                           column(4,numericInput("VGP_ext_cutoff",label = "VGP fixed-filter radius",value = 45))
                                         ),
                                         fluidRow(
                                           column(12,progressBar(
                                             id = "Ext_vgpboot",
                                             value = 0,total=2000,
                                             title = "VGP bootstrap",
                                             display_pct = TRUE))
                                         ),
                                         tableOutput("MVGPpolestat"),
                                         fluidRow(
                                           column(6,actionButton("resetMVGP",label = "Clear current VGP", width = "100%")),
                                           column(6,actionButton("saveMVGP",label = "Add to global list",width = "100%"))
                                         ),
                                         br(),
                                         h4("Preliminary list of loaded VGPs:"),
                                         h5("(Interactive list on Multiple VGPs analysis)"),
                                         br(),
                                         tableOutput("Ext_VGP_list")
                            ),
                            mainPanel(width = 7,
                                      fluidRow(
                                        downloadButton("Ext_VGP_G","Export graph"),
                                        downloadButton("VGP_site_stat","Export current site stat")
                                      ),
                                      column(1),
                                      plotOutput("Ext_VGP_plot")
                            )
                          )
                 ),
                 tabPanel("VGP simulator",
                          sidebarLayout(
                            sidebarPanel(width = 5,
                                         fluidRow(
                                           column(12,h4("VGP (Fisherian) simulator")),
                                         ),
                                         fluidRow(
                                           column(6, textInput("SVGPname",label = "VGP name",value = "VGP_sim")),
                                           column(3,numericInput("SVGPlon",label = "Longitude",value = 0)),
                                           column(3,numericInput("SVGPlat",label = "Latitude",value = 90,max = 90,min = -90))
                                         ),
                                         fluidRow(
                                           column(4,numericInput("SVGPN",label = "N",value = 100)),
                                           column(4,numericInput("SVGPk",label = "K",value = 20)),
                                           column(4,numericInput("k_tol",label = "K tolerance",value = 0.5) %>%
                                                    helper(type = "inline",
                                                           title = "K tolerance",
                                                           content = c("Maximum difference between selected and simulated k. Do not set 0."),
                                                           size = "m",fade = T)),
                                         ),
                                         fluidRow(
                                           column(3,selectInput("SVGP_center",label="Center of plot",
                                                                choices=list("Automatic"= 1,"Manual"= 2), selected=1)),
                                           column(3,numericInput("SVGP_clat",label="Center lat", value=90)),
                                           column(3,numericInput("SVGP_clong",label="Center Long", value=0)),
                                           column(3,selectInput("SVGP_coast", label = "Coastline",
                                                                choices = list("Yes"=1,"no"=2),selected=1)),
                                         ),
                                         fluidRow(
                                           column(3,selectInput("SVGPcolor", label= "VGP color",
                                                                choices= list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,"red"=7,"yellow"=8,"cyan"=9,"gray"=10, "white"=11), selected=2)),
                                           column(3,selectInput("SVGPsymbol", label= "VGP symbol",
                                                                choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1)),
                                           column(3,selectInput("SVGP_stat", label="Statistic",
                                                                choices=list("None"=1,"Fisher"=2,"Bootstrap"=3),selected=2)),
                                           column(3,numericInput("SVGPnb", label = "Bootstrap n.", value = 2000))
                                         ),
                                         fluidRow(
                                           column(12,progressBar(
                                             id = "SVGPboot",
                                             value = 0,total=2000,
                                             title = "VGP bootstrap",
                                             display_pct = TRUE))
                                         ),
                                         tableOutput("SVGPstat"),
                                         fluidRow(
                                           column(6, actionButton("SVGPgo",label = "GENERATE VGP",width = "100%")),
                                           column(6,actionButton("saveSVGP",label = "Add to global list",width = "100%"))
                                         ),
                                         br(),
                                         h4("Preliminary list of loaded VGPs:"),
                                         h5("(Interactive list on Multiple VGPs analysis)"),
                                         br(),
                                         tableOutput("Ext_VGP_list2")

                            ),
                            mainPanel(width = 7,
                                      fluidRow(
                                        downloadButton("SVGP_G","Export graph"),
                                        downloadButton("SVGP_list","Export VGP"),
                                        downloadButton("SVGP_S","Export VGP stat")
                                      ),
                                      column(1),
                                      plotOutput("SVGP_plot")
                            )
                          )
                 ),
                 tabPanel("VGP Rotator",
                          sidebarLayout(
                            sidebarPanel(width = 5,
                                         fluidRow(
                                           column(5,h4("VGPs Euler Rotation")%>%
                                                    helper(type = "inline",
                                                           title = "Info",
                                                           content = c("It works only on one site at a time!!"),
                                                           size = "m",fade = T))
                                         ),
                                         fluidRow(
                                           column(3,selectInput("eul_center",label="Center of plot",
                                                                choices=list("Automatic"= 1,"Manual"= 2), selected=1)),
                                           column(3,numericInput("eul_clat",label="Center lat", value=90)),
                                           column(3,numericInput("eul_clong",label="Center Long", value=0)),
                                           column(3,selectInput("eul_coast", label = "Coastline",
                                                                choices = list("Yes"=1,"no"=2),selected=1)),
                                         ),
                                         fluidRow(
                                           column(3, numericInput("eul_long",label = "EPole Long",value = 0,min = 0,max = 360)),
                                           column(3, numericInput("eul_lat",label = "EPole Lat",value = 90,min = -90,max = 90)),
                                           column(3, numericInput("eul_rot",label = "Rotation",value = 0,min = 0,max = 360)),
                                           column(3, textInput("eul_name",label = "Name",value = "Rotated"))
                                         ),
                                         fluidRow(
                                           column(4,selectInput("eulPlotType",label = "VGPs type",
                                                                choices = list("VGPs"=1,"Fisher"=2,"Bootstrapped"=3),selected = 1)),
                                           column(4,selectInput("eulcolor", label= "VGP color",
                                                                choices= list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,"red"=7,"yellow"=8,"cyan"=9,"gray"=10, "white"=11), selected=8)),
                                           column(4,selectInput("eulsymbol", label= "VGP symbol",
                                                                choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1))
                                         ),
                                         fluidRow(
                                           column(4, actionButton("eulerrot",label = "ROTATE",width = "100%")),
                                           column(4,actionButton("eulerdel",label = "DELETE",width = "100%")),
                                           column(4, actionButton("eulersave",label = "Add to G. List", width = "100%"))
                                         ),
                                         br(),
                                         tableOutput("eul_temp_table"),
                                         fluidRow(
                                           column(5,h4("List of loaded VGPs")%>%
                                                    helper(type = "inline",
                                                           title = "Info",
                                                           content = c("It works only on one site at a time!!"),
                                                           size = "m",fade = T))
                                         ),
                                         br(),
                                         fluidRow(
                                           column(12,DT::dataTableOutput("MVGPlist2"))
                                         ),
                            ),
                            mainPanel(width=7,
                                      fluidRow(
                                        downloadButton("euler_G","Export graph"),
                                        downloadButton("euler_vgp","Export rotated VGP")
                                      ),
                                      column(1),
                                      plotOutput("eulerplot")
                            )
                          )

                 ),
                 tabPanel("Analysis - Page 1 - Average pole",
                          sidebarLayout(
                            sidebarPanel(width = 5,
                                         fluidRow(
                                           column(8,h4("Average pole analysis"))
                                         ),
                                         fluidRow(
                                           column(6,textInput("fileN_MVGP",label = "Average Pole Name",value = "MeanPole"))
                                         ),
                                         fluidRow(
                                           column(3,selectInput("MVGP_center",label="Center of plot",
                                                                choices=list("Automatic"= 1,"Manual"= 2), selected=1)),
                                           column(3,numericInput("MVGP_clat",label="Center lat", value=90)),
                                           column(3,numericInput("MVGP_clong",label="Center Long", value=0)),
                                           column(3,selectInput("MVGP_coast", label = "Coastline",
                                                                choices = list("Yes"=1,"no"=2),selected=1)),
                                         ),
                                         fluidRow(
                                           column(4,selectInput("MVGPsPlotType",label = "VGPs type",
                                                                choices = list("VGPs"=1,"Fisher"=2,"Bootstrapped"=3),selected = 1)),
                                           column(4,selectInput("MVGP_Pole_Stat", label = "M-VGP average",
                                                                choices = list("None"=1,"Fisher on fishers"=2,"Fisher of VGPs"=3,"Bootstrap of VGPs"=4), selected = 1)),
                                           column(4,selectInput("MVGP_names_YN",label = "Plot poles name",
                                                                choices = list("No"=1,"Yes"=2),selected = 1))
                                         ),
                                         fluidRow(
                                           column(4,numericInput("MVGPnb", label = "Bootstrap n.", value = 2000)),
                                           column(4,selectInput("MVGP_aver_sym", label = "Mean sym",
                                                                choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1)),
                                           column(4,selectInput("MVGP_aver_color", label = "Mean color",
                                                                choices= list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,"red"=7,"yellow"=8,"cyan"=9,"gray"=10,"white"=11), selected=10))
                                         ),
                                         fluidRow(
                                           column(12,progressBar(
                                             id = "Mvgpboot",
                                             value = 0,total=2000,
                                             title = "VGP bootstrap",
                                             display_pct = TRUE))
                                         ),
                                         tableOutput("MVGP_ALLVGPS_stat"),
                                         br(),
                                         fluidRow(
                                           column(6,h4("List of loaded VGPs")),
                                           column(6,actionButton("deletevgp",label = "DELETE SELECTED VGP",width = "100%"))
                                         ),
                                         br(),
                                         fluidRow(
                                           column(12,DT::dataTableOutput("MVGPlist"))
                                         ),
                            ),
                            mainPanel(width=7,
                                      fluidRow(
                                        downloadButton("MVGP_G","Export graph"),
                                        downloadButton("VGP_ALL_stat","Export selected VGP average"),
                                        downloadButton("MVGP_list",label = "Export Pole list")
                                      ),
                                      column(1),
                                      plotOutput("MVGP_plot")
                            )
                          )
                 ),
                 tabPanel("Analysis - Page 2 - add paleomag. poles",
                          sidebarLayout(
                            sidebarPanel(width = 5,
                                         fluidRow(
                                           column(8,h4("External poles list"))
                                         ),
                                         fluidRow(
                                           column(6,fileInput("extrapolelist",label = "Load pole list")),
                                           column(6,selectInput("addextreanames",label = "Plot pole names",
                                                                choices = list("No"=1,"Yes"=2),selected = 1))
                                         ),
                                         fluidRow(
                                           column(8,h4("Pole manual entry"))
                                         ),
                                         fluidRow(
                                           column(4,numericInput("extrapolelong",label = "Pole long",value = NULL,min = -360,max = 360)),
                                           column(4,numericInput("extrapolelat",label = "Pole lat",value = NULL,min = -90,max = 90)),
                                           column(4,numericInput("extrapoleA95",label = "95% confidence",value = NULL,min = 0,max = 180))
                                         ),
                                         fluidRow(
                                           column(4,selectInput("extrapolecolor", label= "Pole color",
                                                                choices= list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,"red"=7,"yellow"=8,"cyan"=9,"gray"=10, "white"=11), selected=2)),
                                           column(4,selectInput("extrapolesimbol", label= "Pole symbol",
                                                                choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1)),
                                           column(4,textInput("extrapolename",label = " Manual pole name",value = "E_pole"))
                                         ),
                                         fluidRow(
                                           column(6, actionButton("addextrapole",label = "ADD TO EXTERNAL POLES LIST",width = "100%")),
                                           column(6, actionButton("delextrapolelist",label = "IGNORE LOADED FILE",width = "100%"))
                                         ),
                                         br(),
                                         fluidRow(
                                           column(8,h4("Fisher average of all selected poles"))
                                         ),
                                         fluidRow(
                                           column(3,selectInput("extrapolesfisher",label = "Statistic",
                                                                choices = list("None"=1,"Fisher"=2),selected = 1)),
                                           column(3,selectInput("extrameansymbol",label = "Average symbol",
                                                                choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected = 1)),
                                           column(3,selectInput("extrameancolor", label= "Average color",
                                                                choices= list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,"red"=7,"yellow"=8,"cyan"=9,"gray"=10, "white"=11), selected=2)),
                                           column(3,textInput("extreameanname",label = "Name",value = "Ext_F_mean")),
                                         ),
                                         fluidRow(
                                           column(8,tableOutput("extpolesfisher")),
                                           column(4,actionButton(inputId = "add_F_2_ext_list",label = "Add Fisher to list",width = "100%"))
                                         ),
                                         br(),
                                         fluidRow(
                                           column(6,h4("List of external poles")),
                                           column(6, actionButton("deleteextrapole",label = "DELETE SELECTED POLES",width = "100%"))
                                         ),
                                         br(),
                                         fluidRow(
                                           column(12,DT::dataTableOutput("EP_list"))
                                         ),
                            ),
                            mainPanel(width = 7,
                                      fluidRow(
                                        downloadButton("MVGP_G2",label = "Export graph"),
                                        downloadButton("Ext_pole_fisher_S",label = "Export Fisher stat")
                                      ),
                                      column(1),
                                      plotOutput("MVGP_plot2"))
                          )
                 ),
                 tabPanel("Analysis - Page 3 - APWP",
                          sidebarLayout(
                            sidebarPanel(width = 5,
                                         fluidRow(
                                           column(8,h4("Add APWP"))
                                         ),
                                         fluidRow(
                                           column(4,selectInput("APWP", label = "APWP",
                                                                choices = list("None"=1,"V2023"=2,"T2012"=3,"Custom"=4),selected=1)),
                                           column(4,selectInput("frameV23", label= "V2023 frames",
                                                                choices= list("South Africa"=1,"North America"=2,
                                                                              "South America"=3,"Europe"=4,
                                                                              "India"=5,"Australia"=6,"Antarctica"=7,
                                                                              "Pacific (0-80Ma)"=8,"Iberia (0-80Ma)"=9), selected=1)),
                                           column(4,selectInput("frameT12", label= "T2012 frames",
                                                                choices= list("South Africa"= 1,"North America"=2,
                                                                              "Europe"=3,"India"=4,"Amazonia"=5,
                                                                              "Australia"=6,"East Antarctica"=7), selected=1))
                                         ),
                                         fluidRow(
                                           column(4,fileInput("customAPWP",label = "Custom APWP")),
                                           column(4,numericInput("apwp_Y",label = "APWP min age",value = 0)),
                                           column(4,numericInput("apwp_O",label = "APWP max age",value = 320))
                                         )
                            ),
                            mainPanel(
                              width = 7,
                              fluidRow(
                                downloadButton("MVGP_G3",label = "Export graph")
                              ),
                              column(1),
                              plotOutput("MVGP_plot3")
                            )
                          )
                 ),
               )
      ),
      tabPanel("Magnetic polarity",
               sidebarLayout(
                 sidebarPanel(width = 3,
                              fluidRow(
                                column(12,textInput("fileN_mgstr",label = "Export name",value = "Site"))
                              ),
                              fluidRow(
                                column(6,numericInput("baseMS",label = "Base plot",value = NULL)),
                                column(6,numericInput("topMS",label = "Top plot",value = NULL))
                              ),
                              fluidRow(
                                column(6,numericInput("hGrid",label = "Vert. grid (m)",value = 10)),
                                column(6,selectInput("colmgstr", label= "Point color",
                                                     choices = list("black"=1,"blue"=2,"green"=3,"pink"=4,"purple"=5,"brown"=6,
                                                                    "red"=7,"yellow"=8,"cyan"=9,"gray"=10,"white"=11),selected = 7))
                              ),
                              fluidRow(
                                column(6,numericInput("Doffset",label = "Decl. offset (°)",value = 0))
                              ),
                              br(),
                              fluidRow(
                                column(6,downloadButton("mgstr",label= "Export graph")),
                                column(6,downloadButton("revTab",label="Export Rev."))
                              )
                 ),
                 mainPanel(
                   plotOutput("magstrat")
                 )
               )
      ),
      tabPanel("Site map",
               sidebarLayout(
                 sidebarPanel(width = 3,
                              fluidRow(
                                column(12,textInput("siteText", label= "Site name",value = ""))),
                              fluidRow(
                                column(6,selectInput("siteCol", label= "Color",
                                                     choices = list("black"=1, "blue" = 2, "red" = 3,"green"=4,"purple"=5,"white"=6),selected = 3)),
                                column(6,selectInput("siteSym", label= "Symbol",
                                                     choices = list("circle"=1, "square"=2, "diamond"=3,"triangle"=4),selected=1))),
                              br(),
                              fluidRow(
                                column(12,fileInput("sitefile", label= "Sites file"))
                              ),
                              fluidRow(

                                column(6,selectInput("gridSpace",label="Map grid",
                                                     choices=list("No grid"=1,"10°"=2,"15°"=3,"30°"=4,
                                                                  "45°"=5, "90°"=6), selected=4)),
                                column(6,selectInput("gridCent", label= "Center meridian",
                                                     choices = list("Greenwich"=1, "Anti Greenwich"=2),selected=1))

                              ),
                              fluidRow(
                                column(6,selectInput("landCol", label= "Land color",
                                                     choices = list("black"=1, "gray" = 2,"light gray"=3, "green" = 4,"darkgreen"=5,"light brown"=6,"brown"=7),selected = 3)),
                                column(6,selectInput("seaCol", label= "Sea Color",
                                                     choices = list("cyan"=1, "light cyan"=2, "light green"=3,"white"=4,"light gray"=5),selected=4))
                              ),
                              fluidRow(
                                column(6,selectInput("gridCol", label= "Grid color",
                                                     choices = list("black"=1, "gray" = 2, "light gray" = 3,"blue"=4,"light blue"=5),selected = 2))
                              ),
                              br(),
                              fluidRow(
                                column(12,actionButton("mapgo", label= "Plot - Refresh",width = "100%"))
                              ),
                              br(),
                              fluidRow(
                                column(12,actionButton("resetsitesfile",label = "Delete sites file", width = "100%"))
                              )
                 ),
                 mainPanel(
                   fluidRow(
                     column(3,downloadButton("mapG","Export map"))
                   ),
                   fluidRow(
                     column(12,plotOutput("geomap"))
                   )
                 )
               )
      )
    )
  )



  server <- function(input, output){

    #activate helpers
    observe_helpers()

    ############ VECTOR END-POINTS MODULE

    #read file if present, reset if requested
    #reads three formats for now: Web-DiR, Lamont, Bremen (.cor, one coordinate system)
    sample_list <- reactive({
      if (is.null(input$All_Zijd)){
        return(NULL)
      }
      if(input$Zijd_f_type==1){
        read.csv(file = input$All_Zijd$datapath)
      }
      else if(input$Zijd_f_type==2){
        dat <- read.table(input$All_Zijd$datapath,header = F,skip = 6)
        dat <- dat[,-c(3,4)]
        dat_PmagDiR <- data.frame(matrix(ncol=11,nrow = nrow(dat)))
        colnames(dat_PmagDiR) <- c("sample","step","Sx","Sy","Sz","Gx","Gy","Gz","Bx","By","Bz")
        dat_PmagDiR[,1:2] <- dat[,1:2]
        dat_PmagDiR[,3:5] <- PmagDiR::s2c(DI = dat[,4:5],J = dat[,3])
        dat_PmagDiR[,6:8] <- PmagDiR::s2c(DI = dat[,6:7],J = dat[,3])
        dat_PmagDiR[,9:11] <- PmagDiR::s2c(DI = dat[,8:9],J = dat[,3])
        return(dat_PmagDiR)
      }else if(input$Zijd_f_type==3){
        dat <- read.table(input$All_Zijd$datapath,header = F,skip = 1)
        dat_PmagDiR <- data.frame(matrix(ncol = 11,nrow = nrow(dat)))
        colnames(dat_PmagDiR) <- c("Sample","Step","Sx","Sy","Sz","Gx","Gy","Gz","Bx","By","Bz")
        dat_PmagDiR[,1:2] <- dat[,1:2]
        dat_PmagDiR[,3:5] <- dat[,3:5]
        dat_PmagDiR[,6:8] <- dat[,3:5]
        dat_PmagDiR[,9:11] <- dat[,3:5]
        #convert to A/m
        dat_PmagDiR[,3:11] <- (dat_PmagDiR[,3:11])/0.008
        dat_PmagDiR <- unique(dat_PmagDiR)
        return(dat_PmagDiR)
      }else if(input$Zijd_f_type==4){
        dat <- read.csv(input$All_Zijd$datapath)
        temp_file <- data.frame(matrix(ncol=12,nrow = nrow(dat)))
        colnames(temp_file) <- c("CSF_A","sample","step","Sx","Sy","Sz","Gx","Gy","Gz","Bx","By","Bz")
        temp_file[,1] <- dat[,9]
        temp_file[,2] <- paste(dat[,1],paste(dat[,2],dat[,3],sep=""),
                               paste(dat[,4],dat[,5],sep = ""),
                               paste(dat[,6],dat[,7],sep = ""),
                               paste(dat[,8],(dat[,8]+2),sep = "/"),sep="-")
        temp_file[,3] <- dat[,20]
        temp_file[,4:6] <- PmagDiR::s2c(DI = dat[,12:11],J = dat[,15])
        temp_file[,7:9] <- PmagDiR::s2c(DI = dat[,14:13],J = dat[,15])
        temp_file[,10:12] <- PmagDiR::s2c(DI = dat[,14:13],J = dat[,15])
        temp_file <- temp_file[order(temp_file[,3]),]
        temp_file <- temp_file[order(temp_file[,1]),]
        dat_PmagDiR <- temp_file[,-1]
        return(dat_PmagDiR)
      }else if(input$Zijd_f_type==5){
        dat_PmagDiR <- data.frame(matrix(ncol = 11,nrow = 0))
        colnames(dat_PmagDiR) <- c("Sample","Step","Sx","Sy","Sz","Gx","Gy","Gz","Bx","By","Bz")
        for(i in 1:length(input$All_Zijd[,1])){
          #read first and row and count columns
          f_row <- read.table(input$All_Zijd[[i, 'datapath']],header = F,skip = 2,nrows = 1)
          s_row <- read.table(input$All_Zijd[[i, 'datapath']],header = F,skip = 3,nrows = 1)
          if(ncol(f_row)<ncol(s_row)){
            f_row <- cbind(f_row[1,1],NA,f_row[1,2:ncol(f_row)])
            dat <- read.table(input$All_Zijd[[i, 'datapath']],header = F,skip = 3)
            colnames(f_row) <- colnames(dat)
            dat <- rbind(f_row,dat)
          }else if(ncol(f_row)==ncol(s_row)){dat <- read.table(input$All_Zijd[[i, 'datapath']],header = F,skip = 2)}
          specimen <- input$All_Zijd[[i, 'name']]
          temp_file <- data.frame(matrix(ncol = 11,nrow = nrow(dat)))
          colnames(temp_file) <- c("Sample","Step","Sx","Sy","Sz","Gx","Gy","Gz","Bx","By","Bz")
          temp_file[,1] <- rep(specimen)
          temp_file[,2] <- dat[,2]
          temp_file[,3:5] <- PmagDiR::s2c(DI = dat[,9:10],J = dat[,7])
          temp_file[,6:8] <- PmagDiR::s2c(DI = dat[,3:4],J = dat[,7])
          temp_file[,9:11] <- PmagDiR::s2c(DI = dat[,5:6],J = dat[,7])
          dat_PmagDiR <- rbind(dat_PmagDiR,temp_file)
        }
        return(dat_PmagDiR)
      }
    })

    #create reactive file
    specim <- reactiveValues(specim=NULL)

    #isolate specimen names and make table
    specimens <- reactive({
      if(is.null(sample_list())==F){
        dat <- sample_list()
        specimens <- data.frame(unique(dat[,1]))
        colnames(specimens) <- "specimens"
        specim$list <- specimens
      }
    })

    #send table to UI
    output$samples_list <- DT::renderDataTable(specimens(), server = F,
                                               selection="single",
                                               rownames=F,
                                               options=list(searching=F))

    #isolate specimen selected on list side of the Zijd.
    isolated_specimen <- reactive({
      #depopulate stat temp result files
      specim$DiR_f <- NULL
      specim$DiR_p <- NULL
      specim$DiR_da <- NULL
      specim$DiR_df <- NULL
      specim$DiR_doi <- NULL
      samp <- input$samples_list_rows_selected
      if(length(samp)){
        req(sample_list())
        dat <- sample_list()
        specim$specim <- dat[dat[,1]==specim$list[samp,1],]
      }
    })

    #depopulate stat temp result files when new sample is selected
    observeEvent(input$samples_list_rows_selected,{
      #depopulate stat temp result files
      specim$DiR_f <- NULL
      specim$DiR_p <- NULL
      specim$DiR_da <- NULL
      specim$DiR_df <- NULL
      specim$DiR_doi <- NULL
    })

    #restore specimen
    observeEvent(input$restore_VEPs,{
      specim$specim <- isolated_specimen()
      specim$selectedVEP <- NULL
      specim$selectedVEP <- NULL
      specim$selectedVEP_t <- NULL
      specim$selectedVEP_BW <- NULL
      specim$selectedVEP_BB <- NULL
      specim$selectedVEP_stereo <- NULL
      specim$DiR_f <- NULL
      specim$DiR_p <- NULL
      specim$DiR_da <- NULL
      specim$DiR_df <- NULL
      specim$DiR_doi <- NULL
      specim$DiR_d <- NULL
    })

    #table data of VEP
    output$sampledat <- DT::renderDataTable({
      #if no samples are selected returns null
      samp <- input$samples_list_rows_selected
      if(length(samp)){
        dat <- data.frame(specim$specim[,2])
        colnames(dat) <- "Step"
        dat
      }else{NULL}
    }, server = FALSE, rownames=FALSE,options=list(dom='t',sort=FALSE,
                                                   "drawCallback" = JS("function(settings) {var table = this.api();table.rows().nodes().to$().css('font-size', '12px');}"),
                                                   paging=FALSE),class=list(stripe=FALSE))


    #selection of vector end-points
    selectedVEP <- reactive({
      #select points from step list
      specim$selectedVEP_t <- input$sampledat_rows_selected
      if(input$Zijd_Stereo_shift==1 || input$Zijd_Stereo_shift==2){
        specim$selectedVEP_stereo <- NULL
        #brush white points
        selectedVEP_BW_temp <- rownames(brushedPoints(specim$specim_no_OM, input$plot_brush, xvar = "x", yvar = "z"))
        specim$selectedVEP_BW <- which(rownames(specim$specim) %in% selectedVEP_BW_temp)
        #brush black points
        selectedVEP_BB_temp <- rownames(brushedPoints(specim$specim_no_OM, input$plot_brush, xvar = "x", yvar = "y"))
        specim$selectedVEP_BB <- which(rownames(specim$specim) %in% selectedVEP_BB_temp)
      } else if(input$Zijd_Stereo_shift==3){
        specim$selectedVEP_BW <- NULL
        specim$selectedVEP_BB <- NULL
        #brush stereo points
        selectedVEP_stereo_temp <- rownames(brushedPoints(specim$cart, input$plot_brush, xvar = "x", yvar = "y"))
        specim$selectedVEP_stereo <- which(rownames(specim$specim) %in% selectedVEP_stereo_temp)
      }
      #create single file with all selected points
      specim$selectedVEP <- sort(unique(c(specim$selectedVEP_t,specim$selectedVEP_BW,specim$selectedVEP_BB,specim$selectedVEP_stereo)))
      return(specim$selectedVEP)
    })

    #delete steps
    observeEvent(input$del_VEPs,{
      todelete <- selectedVEP()
      specim$specim <- specim$specim[-todelete,]
      specim$selectedVEP <- NULL
      specim$selectedVEP_t <- NULL
      specim$selectedVEP_BW <- NULL
      specim$selectedVEP_BB <- NULL
      specim$selectedVEP_stereo <- NULL
    })

    #define size of plot, or equal area is too big
    size_plot <- reactive({
      if(input$Zijd_Stereo_shift==1 || input$Zijd_Stereo_shift==2){size <- 800}
      if(input$Zijd_Stereo_shift==3){size <- 600}
      return(size)
    })

    #send Vector end point or equal area fig to UI
    output$zijderveld <- renderPlot({
      #does not send fig if file is not selected
      req(isolated_specimen())
      coord <- input$VEPcoordinates


      #plot Zijderveld
      if(input$Zijd_Stereo_shift==1 || input$Zijd_Stereo_shift==2){
        if(input$VEPticks==1){d_tick=0.05}
        if(input$VEPticks==2){d_tick=0.1}
        if(input$VEPticks==3){d_tick=0.25}
        if(input$VEPticks==4){d_tick=0.5}
        if(input$VEPticks==5){d_tick=1.0}
        ticks <- TRUE
        if(input$VEPticks==6){ticks=FALSE}

        #function plotting the vector end point diagram and return values without order of magnitude
        zijderveld <- function(specim,selected_steps,coordinates=1,orient=1,d_tick=0.25,ticks=T){

          if(coordinates==1){specim <- specim[,c(1,2,3,4,5)]}
          else if(coordinates==2){specim <- specim[,c(1,2,6,7,8)]}
          else if(coordinates==3){specim <- specim[,c(1,2,9,10,11)]}

          colnames(specim) <- c("sample","step","x","y","z")
          #invert z for plotting correct positive and negative
          specim[,5] <- -specim[,5]
          #fix data if North is right and W is up
          if(orient==1){
            specim[,4] <- -specim[,4]
          }
          #fix data if North is up and E is right
          if(orient==2){
            temp <- specim[,3]
            specim[,3] <- specim[,4]
            specim[,4] <- temp
          }
          #eliminate order of magnitudes from values
          OM <- log10(max(abs(specim[,3:5])))
          OM <- ifelse(OM<0,floor(OM), ceiling(OM))
          specim[,3:5] <- (specim[,3:5])/(10**OM)
          #save order of magnitude in reactive file

          #define min and max of x
          if(max(specim[,3])<0 && min(specim[,3])<0){
            x_max <- 0
            x_min <-  min(specim[,3])
            #x_min <- floor(min(specim[,3]))
          }else if(max(specim[,3])>0 && min(specim[,3])>0){
            x_max <-  max(specim[,3])
            x_min <- 0
          }else{
            x_max <- max(specim[,3])
            x_min <- min(specim[,3])
          }
          #define min and max of y
          if(max(specim[,4:5])<0 && min(specim[,4:5])<0){
            y_max <- 0
            #y_min <- floor(min(specim[,4:5]))
            y_min <- min(specim[,4:5])
          }else if(max(specim[,4:5])>0 && min(specim[,4:5])>0){
            y_max <-  max(specim[,4:5])
            y_min <- 0
          }else{
            y_max <- max(specim[,4:5])
            y_min <- min(specim[,4:5])
          }

          #plot empty diagram
          plot(NA,asp=1,
               xlim=c(x_min,x_max),
               ylim=c(y_min,y_max),
               axes=F,xlab="",ylab="",
               #add axes
               panel.first= c(arrows(x0 = x_min,y0 = 0,
                                     x1 = x_max,y1 = 0,length = 0),
                              arrows(x0 = 0,y0 = y_min,
                                     x1 = 0,y1 = y_max, length = 0))
          )
          #ask if ticks are wanted
          if(ticks==T){
            #custom function drawing only the ticks
            draw_ticks <- function(x0, y0, x1, y1){
              t = atan2(y1-y0, x1-x0)
              a1 = pi+t+pi/2
              a2 = pi+t-pi/2
              e1x = x1 + ((x_max-x_min)/100)*cos(a1)
              e1y = y1 + ((x_max-x_min)/100)*sin(a1)
              e2x = x1 + ((x_max-x_min)/100)*cos(a2)
              e2y = y1 + ((x_max-x_min)/100)*sin(a2)
              lines(c(e1x,x1,e2x),c(e1y,y1,e2y))
            }
            #define number of ticks on x and y
            #d_tick is the interval of each unit (*10^OM) subdivision
            x_ticks_l <- floor((-x_min)/d_tick)
            x_ticks_r <- floor((x_max)/d_tick)
            for(i in 0:x_ticks_l){
              draw_ticks(x0 = (-((i-1)*d_tick)),y0 = 0,
                         x1 = (-((i)*d_tick)),y1 = 0)
            }
            for(i in 0:x_ticks_r){
              draw_ticks(x0 = (((i-1)*d_tick)),y0 = 0,
                         x1 = (((i)*d_tick)),y1 = 0)
            }

            y_ticks_d <- floor((-y_min)/d_tick)
            y_ticks_u <- floor((y_max)/d_tick)
            for(i in 0:y_ticks_d){
              draw_ticks(y0 = (-((i-1)*d_tick)),x0 = 0,
                         y1 = (-((i)*d_tick)),x1 = 0)
            }
            for(i in 0:y_ticks_u){
              draw_ticks(y0 = (((i-1)*d_tick)),x0 = 0,
                         y1 = (((i)*d_tick)),x1 = 0)
            }
          }

          #plots vector end points
          points(x=specim[,3],y=specim[,4],type="o",pch=21,bg="black",cex=1.2)
          points(x=specim[,3],y=specim[,5],type="o",pch=21,bg="white",cex=1.2)
          #plot NRM as square
          points(x=specim[1,3],y=specim[1,4],pch=22,bg="black",cex=1.5)
          points(x=specim[1,3],y=specim[1,5],pch=22,bg="white",cex=1.5)
          #highlight selected steps
          if(length(selected_steps)){
            Ssteps <- specim[selected_steps,]
            points(x=Ssteps[,3],y=Ssteps[,4],pch=21,bg="red",cex=1.5)
            points(x=Ssteps[,3],y=Ssteps[,5],pch=21,col="red" ,bg="yellow",cex=1.5)
            if(selected_steps[1]==1){
              points(x=specim[1,3],y=specim[1,4],pch=22,bg="red",cex=1.9)
              points(x=specim[1,3],y=specim[1,5],pch=22,col="red",bg="yellow",cex=1.9)
            }
          }


          #add coordinates
          if(orient==1){
            text(x = x_max,y=0,"N", pos=4, cex=1.2)
            text(x = 0,y=y_max,"W, Up", pos=3, cex=1.2)
          }
          if(orient==2){
            text(x = x_max,y=0,"E", pos=4, cex=1.2)
            text(x = 0,y=y_max,"N, Up", pos=3, cex=1.2)
          }
          #create list for returning both data without order of magnitude AND order of magnitude
          result <- list(specim)
          result$OM <- OM
          return(result)
        }


        #save value with no order of magnitude and plot Zijderveld
        Zijdervel_res <- zijderveld(specim = specim$specim,
                                    selected_steps = selectedVEP(),
                                    coordinates = input$VEPcoordinates,
                                    orient = input$Zijd_Stereo_shift,d_tick = d_tick,ticks = ticks)

        #save data without OM
        specim$specim_no_OM <- Zijdervel_res[[1]]
        #save order of magnitude
        specim$OM <- Zijdervel_res[[2]]

        #save coordinates as number otherwise is a character and uses it for interpolation line
        coordinates <- as.numeric(input$VEPcoordinates)

        if(length(specim$DiR_da)||length(specim$DiR_df)||length(specim$DiR_doi)){
          #add interpolation lines
          if(input$anchor==2){
            if(length(specim$DiR_da)){
              if(input$Zijd_Stereo_shift==1){
                m_xy <- -specim$DiR_da[coordinates,8]
                m_xz <- -specim$DiR_da[coordinates,9]
              }else if(input$Zijd_Stereo_shift==2){
                m_xy <- 1/(specim$DiR_da[coordinates,8])
                m_xz <- -specim$DiR_da[coordinates,10]
              }
              curve((m_xy*x),add = T,col="blue", lty=2,lwd=2)
              curve((m_xz*x),add = T,col="blue", lty=2,lwd=2)
            }
          }else if(input$anchor==1){
            if(length(specim$DiR_df)){
              if(input$Zijd_Stereo_shift==1){
                m_xy <- -specim$DiR_df[coordinates,8]
                m_xz <- -specim$DiR_df[coordinates,9]
                #recalculate x0 y0 and z0 removing OM
                x0 <- specim$DiR_df[coordinates,5]/(10**specim$OM)
                y0 <- -specim$DiR_df[coordinates,6]/(10**specim$OM)
                z0 <- -specim$DiR_df[coordinates,7]/(10**specim$OM)
              }else if(input$Zijd_Stereo_shift==2){
                m_xy <- 1/(specim$DiR_df[coordinates,8])
                m_xz <- -specim$DiR_df[coordinates,10]
                #recalculate x0 y0 and z0 removing OM switching axes x and y!!
                x0 <- specim$DiR_df[coordinates,6]/(10**specim$OM)
                y0 <- specim$DiR_df[coordinates,5]/(10**specim$OM)
                z0 <- -specim$DiR_df[coordinates,7]/(10**specim$OM)
              }
              #plot curve corrected for center of mass
              curve(((m_xy*(x-x0))+y0),add = T,col="blue", lty=2,lwd=2)
              curve(((m_xz*(x-x0))+z0),add = T,col="blue", lty=2,lwd=2)
            }
          }else if(input$anchor==3){
            if(length(specim$DiR_doi)){
              if(input$Zijd_Stereo_shift==1){
                m_xy <- -specim$DiR_doi[coordinates,8]
                m_xz <- -specim$DiR_doi[coordinates,9]
                #recalculate x0 y0 and z0 removing OM
                x0 <- specim$DiR_doi[coordinates,5]/(10**specim$OM)
                y0 <- -specim$DiR_doi[coordinates,6]/(10**specim$OM)
                z0 <- -specim$DiR_doi[coordinates,7]/(10**specim$OM)
              }else if(input$Zijd_Stereo_shift==2){
                m_xy <- 1/(specim$DiR_doi[coordinates,8])
                m_xz <- -specim$DiR_doi[coordinates,10]
                #recalculate x0 y0 and z0 removing OM switching axes x and y!!
                x0 <- specim$DiR_doi[coordinates,6]/(10**specim$OM)
                y0 <- specim$DiR_doi[coordinates,5]/(10**specim$OM)
                z0 <- -specim$DiR_doi[coordinates,7]/(10**specim$OM)
              }
              #plot curve corrected for center of mass
              curve(((m_xy*(x-x0))+y0),add = T,col="blue", lty=2,lwd=2)
              curve(((m_xz*(x-x0))+z0),add = T,col="blue", lty=2,lwd=2)
            }else{NULL}
          }
        }
      }
      #work on equal area plot
      else if(input$Zijd_Stereo_shift==3){
        #select coordinates for table to be used by equal area plot
        if(input$VEPcoordinates==3){specim$specim_t <- specim$specim[,c(1,2,9,10,11)]}
        if(input$VEPcoordinates==2){specim$specim_t <- specim$specim[,c(1,2,6,7,8)]}
        if(input$VEPcoordinates==1){specim$specim_t <- specim$specim[,c(1,2,3,4,5)]}

        #functions converting cartesian to spherical
        c2sD <- function(x,y) {(atan2(y,x))*(180/pi)}
        c2sI <- function(x,y,z) {(asin(z/(sqrt((x^2)+(y^2)+(z^2)))))*(180/pi)}
        c2sInt <- function(x,y,z) {sqrt((x^2)+(y^2)+(z^2))}

        #functions converting degree and radians
        d2r <- function(x) {x*(pi/180)}
        r2d <- function(x) {x*(180/pi)}

        #functions converting inc(x) and dec(y) into equal area
        a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
        a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}

        #add spherical coord columns for table
        dat <- specim$specim_t[,-1]
        dat$D <- (round(c2sD(specim$specim_t[,3],specim$specim_t[,4]),digits = 1))%%360
        dat$I <- round(c2sI(specim$specim_t[,3],specim$specim_t[,4],specim$specim_t[,5]),digits = 1)
        dat$Int <- formatC(c2sInt(specim$specim_t[,3],specim$specim_t[,4],specim$specim_t[,5]), format = "e", digits = 2)
        dat <- subset(dat,select = c(1,5,6,7))
        #add x and y coordinates for equal area plot
        dat$x <- a2cx(abs(dat$I),dat$D)
        dat$y <- a2cy(abs(dat$I),dat$D)
        #copy cartesian coordinates to reactive file for equal area use
        specim$cart <- dat

        #draw base equal area from PmagDiR
        PmagDiR::equalarea()

        #graphical function connecting two points on a sphere with great circle segment
        connect_GC <- function(DI){
          #degrees to radians and vice versa
          d2r <- function(x) {x*(pi/180)}
          r2d <- function(x) {x*(180/pi)}
          #functions converting inc(x) and dec(y) into equal area
          a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
          a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
          #functions spherical (Dec=x, Inc=y) to Cartesian
          s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
          s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
          s2cz <- function(y) {sin(d2r(y))}
          for(i in 2:nrow(DI)){
            #data are data frame 2X2 with dec and inc
            data <- DI[(i-1):i,]
            colnames(data) <- c("dec", "inc")

            ##NEXT PART CALCULATES POLE OF GREAT CIRCLE
            #directions in Cartesian coordinates
            data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
            data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
            data$z <- sin(d2r(data$inc))
            #averaged Cartesian coordinates
            x_av <- mean(data$x)
            y_av <- mean(data$y)
            z_av <- mean(data$z)
            #elements of the distribution matrix
            T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                            sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                            sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
            #distribution matrix
            T <- matrix(T_elements,nrow=3, byrow=TRUE)
            #calculate and copy eigenvalues and vectors
            T_e <- eigen(T)
            T_vec <- T_e$vectors
            #coordinate of V3
            V3inc <- r2d(asin(T_vec[3,3]/(sqrt((T_vec[1,3]^2)+(T_vec[2,3]^2)+(T_vec[3,3]^2)))))
            V3dec <- (r2d(atan2(T_vec[2,3],T_vec[1,3])))%%360
            if(V3inc<0){
              V3dec <- (V3dec+180)%%360
              V3inc <- abs(V3inc)
            }
            #place both points on the horizontal plane
            data_horiz <- PmagDiR::bed_DI(DI = data[,1:2],in_file = FALSE,bed_az = (V3dec+180)%%360,
                                          bed_plunge = (90-V3inc))
            if(max(data_horiz[,1])-min(data_horiz[,1])<180){
              GCS_h <- data.frame(seq(min(data_horiz[,1]),max(data_horiz[,1]),0.5))
              colnames(GCS_h) <- "dec"
              GCS_h$inc <- rep(0)
            }else if((max(data_horiz[,1])-min(data_horiz[,1])>=180)){
              GCS_h <- data.frame((seq(-(360-max(data_horiz[,1])),min(data_horiz[,1]),0.5)))
              colnames(GCS_h) <- "dec"
              GCS_h$inc <- rep(0)
            }
            #place points and connecting great circle back in position
            GCS <- PmagDiR::bed_DI(DI = GCS_h,in_file = F,bed_az = V3dec,bed_plunge = (90-V3inc))
            colnames(GCS) <- c("dec","inc")
            #separate hemispheres
            GCS$em <- ifelse(GCS$inc<=0,-1,1)
            GCS$inc <- abs(GCS$inc)
            GCS$x <- a2cx(GCS$inc,GCS$dec)
            GCS$y <- a2cy(GCS$inc,GCS$dec)
            GCS_upper <- dplyr::filter_all(GCS,all_vars(GCS$em<0))
            GCS_lower <- dplyr::filter_all(GCS,all_vars(GCS$em>0))
            #plot great circles
            points(x = GCS_upper$x,y=GCS_upper$y,type="l", col="black",lty=2)
            points(x = GCS_lower$x,y=GCS_lower$y,type="l", col="black")
          }
        }

        #first draw connecting cirlces
        connect_GC(specim$cart[,2:3])

        #add points
        points(x=specim$cart[,5],
               y=specim$cart[,6],pch=21, bg=ifelse(specim$cart[,3]<0,"white","black"),
               cex=1.4)
        points(x=specim$cart[1,5],
               y=specim$cart[1,6],pch=22, bg=ifelse(specim$cart[1,3]<0,"white","black"),
               cex=1.8)
        #highlight selected points
        if(length(selectedVEP())){
          Hightlight <- specim$cart[selectedVEP(),]
          points(x=Hightlight[,5],
                 y=Hightlight[,6],pch=21, col=ifelse(Hightlight[,3]<0,"red","black"),
                 bg=ifelse(Hightlight[,3]<0,"yellow","red"),
                 cex=1.6)
          if(selectedVEP()[1]==1){
            points(x=Hightlight[1,5],
                   y=Hightlight[1,6],pch=22, col=ifelse(Hightlight[1,3]<0,"red","black"),
                   bg=ifelse(Hightlight[1,3]<0,"yellow","red"),
                   cex=2)
          }
        }
        #plot fisher stat
        if(length(specim$DiR_f)){
          PmagDiR::plot_a95(specim$DiR_f[coord,1],specim$DiR_f[coord,2],specim$DiR_f[coord,3],
                            on_plot = T,symbol = "d",col_d = "purple",col_u = "pink")
        }
        if(length(specim$DiR_p)){
          specim$DiR_p <- as.data.frame(specim$DiR_p)
          PmagDiR::plot_plane(specim$DiR_p[coord,1],specim$DiR_p[coord,2],on_plot = T,col_cU = "blue",col_cD = "blue",symbol = "d",col_d = "purple")
        }
      }

      #save plot for export figure
      specim$savedplot <- recordPlot()
    }, height = reactive({size_plot()}))

    #creates reactive file for saving steps of PCA
    specim$saved_steps <- data.frame(matrix(ncol=1,nrow=0))

    #calculate PCA or fisher and save steps selection
    observeEvent(input$runVEPstat,{
      req(specim$specim)
      req(length(selectedVEP())>1)

      DiR <- NULL
      c <- input$anchor
      if(c==1 || c==2 || c==3 || c==5){
        #calculate PCA-derived direction and MAD from demagnetization steps
        #VEPs is expressed in Cartesian coordinates x,y,z
        run_PCA <- function(VEPs,anchor) {
          #degree to radians and VV
          d2r <- function(x) {x*(pi/180)}
          r2d <- function(x) {x*(180/pi)}
          data <- VEPs
          colnames(data) <- c("x", "y","z")
          #averaged Cartesian coordinates
          x_av <- mean(data$x)
          y_av <- mean(data$y)
          z_av <- mean(data$z)
          #copy coordinates for anchored directions or great circle
          if(anchor==1) {
            #calculate coordinates with new center of mass for PCA
            data$xn <- data$x-x_av
            data$yn <- data$y-y_av
            data$zn <- data$z-z_av
          }
          else if (anchor==2){
            data$xn <- data$x
            data$yn <- data$y
            data$zn <- data$z
          }
          else if(anchor==3) {
            #includes origin and calculate new center of mass
            newrow <- c(0,0,0)
            data <- rbind(data,newrow)
            data$xn <- data$x-x_av
            data$yn <- data$y-y_av
            data$zn <- data$z-z_av
          }
          else if (anchor==5){
            #if great circle, data must be transformed in unit vectors as suggested by MF ME 1988 (EPSL87)
            #I use the PmagDiR::c2s and s2c funtions, made for dec inc, and eliminating vector length
            data_spherical <- PmagDiR::c2s(data)
            data <- PmagDiR::s2c(data_spherical)
            data$xn <- data$x
            data$yn <- data$y
            data$zn <- data$z
          }
          #elements of the distribution matrix
          T_elements <- c(sum((data$xn)*(data$xn)),sum(data$xn*data$yn),sum(data$xn*data$zn),
                          sum(data$yn*data$xn),sum(data$yn*data$yn),sum(data$yn*data$zn),
                          sum(data$zn*data$xn),sum(data$zn*data$yn),sum(data$zn*data$zn))

          Tm <- matrix(T_elements,3, 3)
          T_e <- eigen(Tm)
          T_vec <- T_e$vectors
          T_val <- T_e$value

          #interpolate line through points
          if(anchor==1 || anchor== 2 || anchor==3){
            #calculate dec inc of max variance
            Vdec <- (r2d(atan2(T_vec[2,1],T_vec[1,1])))%%360
            Vinc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))

            #flipping V1 module, if directions goes opposite to vector tip
            tip <- c(data[1,1]-data[nrow(data),1],data[1,2]-data[nrow(data),2],data[1,3]-data[nrow(data),3])
            tipdec <- (r2d(atan2(tip[2],tip[1])))%%360
            tipinc <- r2d(asin(tip[3]/(sqrt((tip[1]^2)+(tip[2]^2)+(tip[3]^2)))))
            deltadec_tip_V1<- abs(tipdec-Vdec)
            dist_tip_V1 <- r2d(acos((sin(d2r(tipinc))*sin(d2r(Vinc)))+
                                      (cos(d2r(tipinc))*cos(d2r(Vinc))*cos(d2r(deltadec_tip_V1)))))
            if(dist_tip_V1>90){
              Vdec <- (Vdec+180)%%360
              Vinc <- -Vinc
            }
            #calculate max ang dev of line
            MAD <- r2d(atan(sqrt(((T_val[2])+(T_val[3]))/T_val[1])))


            #calculate x y z coordinates of V1
            V1_x <- cos(d2r(Vdec))*cos(d2r(Vinc))
            #next because zijderveld y axis is down pointing
            V1_y <- (sin(d2r(Vdec))*cos(d2r(Vinc)))
            V1_z <- (sin(d2r(Vinc)))

            #calculate inclination of interpolating lines
            m_xy <- V1_y/V1_x
            m_xz <- V1_z/V1_x
            m_yz <- V1_z/V1_y


          }else if(anchor==5){
            #calculate dec inc of min variance that is the pole of the plane or circle
            Vdec <- (r2d(atan2(T_vec[2,3],T_vec[1,3])))%%360
            Vinc <- r2d(asin(T_vec[3,3]/(sqrt((T_vec[1,3]^2)+(T_vec[2,3]^2)+(T_vec[3,3]^2)))))
            #flip pole if negative
            if(Vinc<0){
              Vdec <- (Vdec+180)%%360
              Vinc <- abs(Vinc)
            }
            MAD <- r2d(atan(sqrt((T_val[3]/T_val[2])+(T_val[3]/T_val[1]))))
          }

          #number of data points
          N <- nrow(data)

          dirs <- cbind(Vdec,Vinc,MAD,N)
          colnames(dirs) <- c("Dec", "Inc","MAD","N")

          #add cartesian coordinates of V1 for calculating interpolating lines, only if line, and center of mass in original OM
          if(anchor==1 || anchor==2 || anchor==3){
            dirs <- cbind(dirs,x_av,y_av,z_av,m_xy,m_xz,m_yz)
          }
          return(dirs)
        }

        #create file with selected VEPs from GUI
        VEPs_for_PCA_sp <- specim$specim[selectedVEP(),3:5]
        VEPs_for_PCA_geo <- specim$specim[selectedVEP(),6:8]
        VEPs_for_PCA_tc <- specim$specim[selectedVEP(),9:11]

        #performs PCA
        DiR_sp <- run_PCA(VEPs = VEPs_for_PCA_sp, anchor = input$anchor)
        DiR_geo <- run_PCA(VEPs = VEPs_for_PCA_geo, anchor = input$anchor)
        DiR_tc <- run_PCA(VEPs = VEPs_for_PCA_tc, anchor = input$anchor)

      } else if(c==4){
        #degree to radians and VV
        d2r <- function(x) {x*(pi/180)}
        r2d <- function(x) {x*(180/pi)}
        dat_sp <- specim$specim[selectedVEP(),3:5]
        dat_geo <- specim$specim[selectedVEP(),6:8]
        dat_tc <- specim$specim[selectedVEP(),9:11]

        dat_sp$dec <- (r2d(atan2(dat_sp[,2],dat_sp[,1])))%%360
        dat_sp$inc <- r2d(asin(dat_sp[,3]/(sqrt((dat_sp[,1]^2)+(dat_sp[,2]^2)+(dat_sp[,3]^2)))))

        dat_geo$dec <- (r2d(atan2(dat_geo[,2],dat_geo[,1])))%%360
        dat_geo$inc <- r2d(asin(dat_geo[,3]/(sqrt((dat_geo[,1]^2)+(dat_geo[,2]^2)+(dat_geo[,3]^2)))))

        dat_tc$dec <- (r2d(atan2(dat_tc[,2],dat_tc[,1])))%%360
        dat_tc$inc <- r2d(asin(dat_tc[,3]/(sqrt((dat_tc[,1]^2)+(dat_tc[,2]^2)+(dat_tc[,3]^2)))))

        DiR_sp <- PmagDiR::fisher(dat_sp[,4:5],export = F)
        DiR_geo <- PmagDiR::fisher(dat_geo[,4:5],export = F)
        DiR_tc <- PmagDiR::fisher(dat_tc[,4:5],export = F)
      }
      #compile different files for Zijderveld graph
      if(c==1){
        #copy result on reactive file and empty others
        specim$DiR_df <- rbind(DiR_sp,DiR_geo,DiR_tc)
        specim$DiR_da <- NULL
        specim$DiR_doi <- NULL
        specim$DiR_f <- NULL
        specim$DiR_p <- NULL
      }
      else if(c==2){
        specim$DiR_df <- NULL
        specim$DiR_da <- rbind(DiR_sp,DiR_geo,DiR_tc)
        specim$DiR_doi <- NULL
        specim$DiR_f <- NULL
        specim$DiR_p <- NULL
      }
      else if(c==3){
        specim$DiR_df <- NULL
        specim$DiR_da <- NULL
        specim$DiR_doi <- rbind(DiR_sp,DiR_geo,DiR_tc)
        specim$DiR_f <- NULL
        specim$DiR_p <- NULL
      }
      else if(c==4){
        specim$DiR_f <- rbind(DiR_sp,DiR_geo,DiR_tc)
        specim$DiR_df <- NULL
        specim$DiR_da <- NULL
        specim$DiR_doi <- NULL
        specim$DiR_p <- NULL
      }
      else if(c==5){
        specim$DiR_p <- rbind(DiR_sp,DiR_geo,DiR_tc)
        specim$DiR_df <- NULL
        specim$DiR_da <- NULL
        specim$DiR_doi <- NULL
        specim$DiR_f <- NULL
      }

      #save demagnetization steps in a temporary file, used if save is requested
      specim$saved_steps_temp <- paste(specim$specim[selectedVEP(),2],collapse = ",")
    })

    # creates text with results
    PCA_result <- reactive({
      c <- input$anchor
      coordinates <- as.numeric(input$VEPcoordinates)
      #assign same name to different files
      if(c==1) {if(length(specim$DiR_df)){specim$DiR_d <- specim$DiR_df}else{specim$DiR_d <- NULL}}
      if(c==2) {if(length(specim$DiR_da)) {specim$DiR_d <- specim$DiR_da}else{specim$DiR_d <- NULL}}
      if(c==3) {if(length(specim$DiR_doi)) {specim$DiR_d <- specim$DiR_doi}else{specim$DiR_d <- NULL}}
      if(c==4) {if(length(specim$DiR_f)) {specim$DiR_d <- specim$DiR_f}else{specim$DiR_d <- NULL}}
      if(c==5) {if(length(specim$DiR_p)) {specim$DiR_d <- specim$DiR_p}else{specim$DiR_p <- NULL}}
      if(length(specim$DiR_d)){
        if(c==1 || c==2 || c==3 || c==5){
          PCA_text <- paste("N= ",round(specim$DiR_d[coordinates,4],digits = 0),", ",
                            "Declination= ",round(specim$DiR_d[coordinates,1],digits = 2),", ",
                            "Inclination= ",round(specim$DiR_d[coordinates,2],digits = 2),", ",
                            "M.A.D.= ", round(specim$DiR_d[coordinates,3],digits = 2),sep = "")
        }
        else if(c==4){
          PCA_text <- paste("N= ",round(specim$DiR_d[coordinates,4],digits = 0),", ",
                            "Declination= ",round(specim$DiR_d[coordinates,1],digits = 2),", ",
                            "Inclination= ",round(specim$DiR_d[coordinates,2],digits = 2),", " ,
                            "a95= ",round(specim$DiR_d[coordinates,3],digits = 2),", ",
                            "k=  ",round(specim$DiR_d[coordinates,6],digits = 2), sep="")
        }
      }
      else{PCA_text <- NULL}
      return(PCA_text)
    })

    #send results to figure
    output$PCA_result <- renderText({
      req(PCA_result())
      PCA_result()
    })

    #creates reactive result file
    specim$PCA_result_file <- data.frame(matrix(ncol=14,nrow = 0))

    #save results in reactive file
    observeEvent(input$save_PCA,{
      #DiR_d contains the results no matter the interpolation
      req(specim$DiR_d)
      c <- input$anchor
      temp_result <- data.frame(matrix(ncol=14,nrow = 1))
      colnames(temp_result) <- c("Sample","N","S_Dec","S_Inc","G_Dec","G_Inc","B_Dec","B_Inc","MAD","a95","k","Type","Comp","Steps")
      temp_result[1,1] <- unique(specim$specim[,1])
      temp_result[1,2] <- round(specim$DiR_d[1,4],digits = 0)
      temp_result[1,3] <- round(specim$DiR_d[1,1],digits = 2)
      temp_result[1,4] <- round(specim$DiR_d[1,2],digits = 2)
      temp_result[1,5] <- round(specim$DiR_d[2,1],digits = 2)
      temp_result[1,6] <- round(specim$DiR_d[2,2],digits = 2)
      temp_result[1,7] <- round(specim$DiR_d[3,1],digits = 2)
      temp_result[1,8] <- round(specim$DiR_d[3,2],digits = 2)
      temp_result[1,9] <- ifelse(c==4,"",round(specim$DiR_d[1,3],digits = 2))
      temp_result[1,10] <- ifelse(c==4,round(specim$DiR_d[1,3],digits = 2),"")
      temp_result[1,11] <- ifelse(c==4,round(specim$DiR_d[1,6],digits = 2),"")
      if(c==1) {Type <- "PCA_F"}
      if(c==2) {Type <- "PCA_A"}
      if(c==3) {Type <- "PCA_OI"}
      if(c==4) {Type <- "Fisher"}
      if(c==5) {Type <- "GC"}
      temp_result[1,12] <- Type
      temp_result[1,13] <- input$comp_name
      #add steps of all interpreted samples
      temp_result[1,14] <- specim$saved_steps_temp
      specim$PCA_result_file <- rbind(specim$PCA_result_file,temp_result)
    })

    #save interpretation and steps file as .csv text
    output$export_PCA <- downloadHandler(
      filename = function() {
        paste("Interpolated_directions_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(specim$PCA_result_file, file,row.names = FALSE)
      }
    )

    #reactive name of sample for saving figure
    sampleName <- reactive({
      if(input$Zijd_Stereo_shift==1 || input$Zijd_Stereo_shift==2){samplename <- paste(specim$specim[1,1],"_VEP_")}
      else if(input$Zijd_Stereo_shift==3){samplename <- paste(specim$specim[1,1],"_EA_")}
      return(samplename)
    })

    #export figure VEPs or Equal area
    output$export_Zijd <- downloadHandler(
      filename = function() {
        paste(sampleName(), Sys.Date(), ".pdf", sep="")
      },
      content = function(file) {
        pdf(file, onefile = TRUE,width = 10,height = 10)
        replayPlot(specim$savedplot)
        dev.off()
      }
    )

    #create text file with unit
    zijd_unit <- reactive({
      if(input$Zijd_Stereo_shift==3 || input$VEPticks==6 || length(isolated_specimen())==0){
        z_unit <- NULL
      }else if(input$Zijd_Stereo_shift==1 || input$Zijd_Stereo_shift==2){
        if(input$VEPticks==1) fract <- "0.05 E"
        if(input$VEPticks==2) fract <- "0.1 E"
        if(input$VEPticks==3) fract <- "0.25 E"
        if(input$VEPticks==4) fract <- "0.5 E"
        if(input$VEPticks==5) fract <- "1 E"
        z_unit <- paste("Unit: ",fract,specim$OM," ",input$textunit,sep = "")
      }
      return(z_unit)
    })

    #print unit of Vector end points axes
    output$Zijd_Unit <- renderText({
      req(zijd_unit())
      zijd_unit()
    })

    ##### All saved samples page


    #import external file and attche it to internal result file if exists
    observeEvent(input$import_PCA,{
      specim$tab_result_ext <- read.csv(file = input$import_PCA$datapath)
      if(!is.null(specim$PCA_result_file)) {
        specim$PCA_result_file <- rbind(specim$PCA_result_file,specim$tab_result_ext)
      }
    })


    TAB <- reactive({
      if(nrow(specim$PCA_result_file)==0){return(NULL)}
      else{
        result_table <- specim$PCA_result_file[,-c(2,3,4,11,14)]
        return(result_table)
      }
    })

    #turn table in a interactive table
    output$saved_interpol <- DT::renderDataTable({TAB()}, server = FALSE, rownames=F, extension= 'Scroller',
                                                 options=list(dom='t',sort=T, paging=T,
                                                              deferRender = TRUE,
                                                              scrollY = 600,
                                                              scroller = TRUE,
                                                              "drawCallback" = JS("function(settings) {var table = this.api();table.rows().nodes().to$().css('font-size', '12px');}"),
                                                              initComplete = JS(
                                                                "function(settings, json) {",
                                                                "$(this.api().table().header()).css({'font-size': '85%'});",
                                                                "}")),
                                                 class=list(stripe=FALSE))

    #delete selected directions permanently
    observeEvent(input$del_interpol,{
      to_delete <- input$saved_interpol_rows_selected
      if(length(input$saved_interpol_rows_selected)>0){
        specim$PCA_result_file_BU <- specim$PCA_result_file
        specim$PCA_result_file <- specim$PCA_result_file[-to_delete,]
      }
    })

    #undo delete
    observeEvent(input$undel_interpol,{
      if(nrow(specim$PCA_result_file_BU)>0){
        specim$PCA_result_file <- specim$PCA_result_file_BU
      }
    })

    #export only selected directions
    output$export_interpol <- downloadHandler(
      filename = function() {
        paste(input$sel_interpol_name,"_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        sel <- input$saved_interpol_rows_selected
        data <- specim$PCA_result_file[sel,]
        write.csv(data, file,row.names = FALSE)
      }
    )

    #Qui funziona ma voglio metterlo nella parte grafica... forse
    observeEvent(input$comb_DI_GC,{
      req(TAB())
      dat <- input$saved_interpol_rows_selected
      dirs_full <- specim$PCA_result_file
      dirs_selected <- dirs_full[dat,]
      if(length(dat)>0){
        dirs_temp <- dirs_selected[dirs_selected$Type!="GC",]
        if(input$EAcoordinates==1){dirs <- dirs_temp[,5:6]}
        else if(input$EAcoordinates==2){dirs <- dirs_temp[,7:8]}
        circles_temp <- dirs_selected[dirs_selected$Type=="GC",]
        if(input$EAcoordinates==1){circles <- circles_temp[,5:6]}
        else if(input$EAcoordinates==2){circles <- circles_temp[,7:8]}
      }
      if(nrow(circles)>0){
        #create empty directions file if no directions are selected
        if(nrow(dirs)>0){
          DI <- dirs
        }else{
          DI <- data.frame(matrix(ncol = 2,nrow = 0))
          colnames(DI) <- c("dec","inc")
        }
        #save some details to fill table
        specim$GC_firstCols <- data.frame(circles_temp[,1:2])
        specim$GC_lastCols <- data.frame(circles_temp[,9:14])

        #apply function from PmagDiR
        specim$GC_directions <- PmagDiR::comb_GC_dirs(dirs = DI,poles = circles)
      }
    })

    #erase file if coordinate are changed
    observeEvent(input$EAcoordinates,{specim$GC_directions <- NULL})

    #erase GC file if asked
    observeEvent(input$GC_erase,{specim$GC_directions <- NULL})

    #save file to list
    observeEvent(input$save_GC,{
      req(specim$GC_directions)
      #create empty tab
      temp_table <- data.frame(matrix(ncol=14,nrow = nrow(specim$GC_directions)))
      colnames(temp_table) <- c("Sample","N","S_Dec","S_Inc","G_Dec","G_Inc","B_Dec","B_Inc","MAD","a95","k","Type","Comp","Steps")
      #paste details copied above
      temp_table[1:nrow(temp_table),1:2] <- specim$GC_firstCols
      temp_table[1:nrow(temp_table),9:14] <- specim$GC_lastCols
      if(input$EAcoordinates==2){
        temp_table[1:nrow(temp_table),7:8] <- round(specim$GC_directions, digits = 2)
      } else if(input$EAcoordinates==1){
        temp_table[1:nrow(temp_table),5:6] <- round(specim$GC_directions, digits = 2)
      }
      temp_table[,12] <- rep("Dir")
      specim$PCA_result_file <- rbind(specim$PCA_result_file,temp_table)
      specim$GC_directions <- NULL
    })




    #plot directions and circles taking data from table of results
    output$saved_interpol_EA <- renderPlot({
      req(TAB())
      dat <- input$saved_interpol_rows_selected
      dirs_full <- TAB()
      dirs_selected <- dirs_full[dat,]
      #cut great circles temporarily
      dirs_temp <- dirs_selected[dirs_selected$Type!="GC",]
      if(input$EAcoordinates==1){dirs <- dirs_temp[,2:3]}
      else if(input$EAcoordinates==2){dirs <- dirs_temp[,4:5]}
      PmagDiR::plot_DI(dirs)
      #add great circles
      circles_temp <- dirs_selected[dirs_selected$Type=="GC",]
      if(nrow(circles_temp>=1)){
        if(input$EAcoordinates==1){circles <- circles_temp[,2:3]}
        else if(input$EAcoordinates==2){circles <- circles_temp[,4:5]}
        #plot plane is designed for a single circle so it has to be reiterated
        for(i in 1:nrow(circles)){
          PmagDiR::plot_plane(circles[i,1],circles[i,2],on_plot = TRUE,col_cU = "blue",col_cD = "blue",symbol = "d",col_d = "yellow")
        }
      }
      if(is.null(specim$GC_directions)==FALSE){PmagDiR::plot_DI(DI = specim$GC_directions,on_plot = T,
                                                                col_d = "red",col_u = "pink",symbol = "t")}
      #add GC directions estimate if exist
    },height = 700,width = 700)


    ############ MAPPING MODULE
    #creates reactive value for checking if file is uploaded
    values <- reactiveValues(mapsites = NULL)
    #check for uploaded file
    observeEvent(input$sitefile,{values$mapsites <- "uploaded"})
    #reset upload if requested
    observeEvent(input$resetsitesfile,{values$mapsites <- "reset"})
    #read file if present, reset if requested
    sites_file <- reactive({
      if (is.null(values$mapsites)) {
        return(NULL)
      } else if (values$mapsites == 'uploaded') {
        read.csv(file = input$sitefile$datapath)
      } else if (values$mapsites == 'reset') {
        return(NULL)
      }
    })


    geo_point_plot <- eventReactive(input$mapgo,{
      geo_point <- function(){
        d2r <- function(x) {x*(pi/180)}
        r2d <- function(x) {x*(180/pi)}

        #functions converting long & lat to xy in KavrayskiyVII projection
        c2x <- function(lon,lat) {((3*d2r(lon))/2)*(sqrt((1/3)-((d2r(lat)/pi)^2)))}
        c2y <- function(lat) {d2r(lat)}

        #functions spherical (lon=x, lat=y) to Cartesian
        s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
        s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
        s2cz <- function(y) {sin(d2r(y))}

        #define center meridian
        if(input$gridCent==1) {center <- 0}
        else if(input$gridCent==2){center <- 180}
        #define land color
        if(input$landCol==1){landcolor <- "black"}
        if(input$landCol==2){landcolor <- "gray"}
        if(input$landCol==3){landcolor <- "lightgray"}
        if(input$landCol==4){landcolor <- "green"}
        if(input$landCol==5){landcolor <- "darkgreen"}
        if(input$landCol==6){landcolor <- "orange"}
        if(input$landCol==7){landcolor <- "brown"}
        #define sea color
        if(input$seaCol==1){seacolor <- "cyan"}
        if(input$seaCol==2){seacolor <- "lightcyan"}
        if(input$seaCol==3){seacolor <- "lightgreen"}
        if(input$seaCol==4){seacolor <- "white"}
        if(input$seaCol==5){seacolor <- "lightgray"}
        #define grid color
        if(input$gridCol==1){gridcolor <- "black"}
        if(input$gridCol==2){gridcolor <- "gray"}
        if(input$gridCol==3){gridcolor <- "lightgray"}
        if(input$gridCol==4){gridcolor <- "blue"}
        if(input$gridCol==5){gridcolor <- "lightblue"}

        #Grid selection
        if(input$gridSpace==1) {grid <- 0}
        if(input$gridSpace==2) {grid <- 10}
        if(input$gridSpace==3) {grid <- 15}
        if(input$gridSpace==4) {grid <- 30}
        if(input$gridSpace==5) {grid <- 45}
        if(input$gridSpace==6) {grid <- 90}

        #draw map
        Map_KVII(grid=grid,center=center,seaCol = seacolor,landCol = landcolor,gridCol = gridcolor)

        #plot point
        S_lon <- input$long
        S_lon <- S_lon-center
        S_lon <- ifelse(S_lon>180,
                        S_lon-360,S_lon)
        S_lon <- ifelse(S_lon<(-180),S_lon+360,S_lon)
        S_lat <- input$lat
        #select symbol
        if(input$siteSym==1) pch <- 21
        if(input$siteSym==2) pch <- 22
        if(input$siteSym==3) pch <- 23
        if(input$siteSym==4) pch <- 24

        #select color
        if(input$siteCol==1) col <- "black"
        if(input$siteCol==2) col <- "blue"
        if(input$siteCol==3) col <- "red"
        if(input$siteCol==4) col <- "darkgreen"
        if(input$siteCol==5) col <- "purple"

        site_x <- c2x(S_lon,S_lat)
        site_y <- c2y(S_lat)
        points(x=site_x,y=site_y,pch=pch, col="black",bg=col,cex=1.2)
        if(is.null(input$siteText)==FALSE){
          sitename <- input$siteText
          text(x=site_x, y=site_y,pos=3,substitute(paste(bold(sitename))), cex= 1.5)
        }
        #req(sites_file())
        if(is.null(sites_file())==FALSE){
          sites <- sites_file()
          for(i in 1:nrow(sites)){
            lon <- sites[i,2]
            lon <- lon-center
            lon <- ifelse(lon>180,lon-360,lon)
            lon <- ifelse(lon<(-180),lon+360,lon)
            x <- c2x(lon,sites[i,3])
            y <- c2y(sites[i,3])
            if(sites[i,4]=="c") pch <- 21
            if(sites[i,4]=="s") pch <- 22
            if(sites[i,4]=="d") pch <- 23
            if(sites[i,4]=="t") pch <- 24
            points(x=x,y=y,pch=pch, col="black",bg=sites[i,5],cex=1.2)
            name <- sites[i,1]
            text(x=x, y=y,pos=3,substitute(paste(bold(name))), cex= 1.5)
          }
        }
      }
      geo_point()
      mapPlot <- recordPlot()
      output$mapG <- downloadHandler(
        filename = function() {
          paste("Map_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 15,height = 10)
          replayPlot(mapPlot)
          dev.off()
        }
      )

    })
    output$geomap <- renderPlot({
      geo_point_plot()
    },height = 700)
    ############ END OF MAPPING MODULE

    ############ Takes direction input file and fix it depending on commands
    #creates reactive value for checking if file is uploaded
    values <- reactiveValues(Dirs = NULL)
    #check for uploaded file and save name
    observeEvent(input$file,{values$Dirs <- "uploaded"})
    observeEvent(input$file,{
      file <- input$file
      dirsFileName <- file$name
      assign("dirsFileName", dirsFileName, envir = .GlobalEnv)
    })
    #reset upload if requested
    observeEvent(input$resetDir,{values$Dirs <- "reset"})
    #read file if present, reset if requested
    input_file <- reactive({
      if (is.null(values$Dirs)) {
        return(NULL)
      } else if (values$Dirs == 'uploaded') {
        read.csv(file = input$file$datapath)
      } else if (values$Dirs == 'reset') {
        return(NULL)
      }
    })

    #fix dirs coordinate depending on input, file ncol, different cutoff
    fix_DI <- function(input_file,file=input$filetype,
                       Slat=input$lat, Slong=input$long,
                       coord=input$coord, cutoff=input$cutoff,
                       VGP_fixed=input$VGP_fixed, MinInc=input$MinInc, MaxInc=input$MaxInc){
      req(input_file)
      DIRS <- input_file
      if(file==1){
        if(coord==1 || coord==2 || coord==3){
          DI <- DIRS
        }
      }
      if(file==2){
        if(coord==1 || coord==3){DI <- DIRS}
        if(coord==2){DI <- PmagDiR::bed_DI(DIRS)}
      }
      if(file==3){
        #save extra temporary file for using bed coordinates for cutoff but plot Geo coord
        DIBB <- DIRS
        if(coord==1 || coord==3){DI <- DIRS[,-c(3,4)]}
        if(coord==2){DI <- DIRS[,-c(1,2)]}
      }
      if(file==4){                         #NON VA QUESTA OPZIONE!!!!
        #save extra temporary file for using bed coordinates for cutoff but plot Geo coord
        DIBB <- DIRS
        if(coord==1){DI <- DIRS[,5:6]}
        if(coord==2){DI <- DIRS[,7:8]}
        if(coord==3){DI <- DIRS[,3:4]}
      }
      #apply cutoff & filters
      if(file==2 && coord==1){geo=TRUE}else{geo=FALSE}
      #in case of dirsfile type 3 and geo coordinates takes bedding coordinate to use as filter after cutoff
      if(file==3 && coord==1){
        if(cutoff>=2 && cutoff<=5){
          DI <- DIRS[,3:4]
        }
      }
      #in case of dirsfile type webDiR and geo or specimen coordinates takes bedding coordinate to use as filter after cutoff
      if(file==4){
        if(cutoff>=2 && cutoff<=5){
          if(coord==2 || coord==3){DI <- DIRS[,7:8]}
        }
      }
      if(cutoff==2){DI <- PmagDiR::cut_DI(DI = DI,lat=Slat,long = Slong,geo = geo,Shiny = T)}
      else if(cutoff==3){DI <- PmagDiR::cut_DI(DI = DI,lat=Slat,long = Slong,inc_f = F,geo = geo,Shiny = T)}
      else if(cutoff==4){DI <- PmagDiR::cut_DI(DI = DI,VD=F,cutoff = VGP_fixed ,lat=Slat,long = Slong,geo = geo,Shiny = T)}
      else if(cutoff==5){DI <- PmagDiR::cut_DI(DI = DI,VD=F,cutoff = VGP_fixed ,lat=Slat,long = Slong, inc_f=F,geo = geo,Shiny = T)}
      else if(cutoff==6){DI <- DI[DI[,2]>0,]}
      else if(cutoff==7){DI <- DI[DI[,2]<0,]}
      else if(cutoff==8){
        DI1 <- DI[DI[,2]<MinInc,]
        DI2 <- DI[DI[,2]>MaxInc,]
        DI <- rbind(DI1,DI2)
      }
      #in case of dirsfile type 2 and geo coordinates gives back geo coordinate dirs filtered by rownames after cutoff
      if(file==3 && coord==1){DI <- DIBB[rownames(DI),1:2]}
      if(file==4){
        if(coord==2){DI <- DIBB[rownames(DI),7:8]}
        if(coord==3){DI <- DIBB[rownames(DI),3:4]}
      }
      return(DI)
    }

    ############ EQUAL AREA MODULE
    #modified fisher_plot function
    fisher_plot_S <- function(DI, plot=TRUE, col_d="red",col_u="white",col_l="black",symbol="c") {
      d2r <- function(x) {x*(pi/180)}
      r2d <- function(x) {x*(180/pi)}
      data <- DI
      data <- na.omit(data)
      data <- data[,1:2]
      colnames(data) <- c("dec", "inc")
      #directions in Cartesian coordinates
      data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
      data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
      data$z <- sin(d2r(data$inc))
      #averaged Cartesian coordinates
      x_av <- mean(data$x)
      y_av <- mean(data$y)
      z_av <- mean(data$z)
      #elements of the distribution matrix
      T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                      sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                      sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
      #distribution matrix
      T <- matrix(T_elements,nrow=3, byrow=TRUE)
      #calculate and copy eigenvalues and vectors
      T_e <- eigen(T,symmetric = TRUE)
      T_vec <- T_e$vectors
      T_val <- T_e$values
      #calculate dec inc of max variance
      V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
      V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
      V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)
      #next  calculates difference between dec_inc and average
      data$Dec_aver <- rep(V1dec)
      data$Inc_aver <- rep(V1inc)
      data$delta <- abs(data$dec-data$Dec_aver)
      data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                              (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
      #Isolate modes
      if(any(data$diff<=90)){
        mode1 <- as.data.frame(data$dec[data$diff<=90])
        mode1$inc <- data$inc[data$diff<=90]
        colnames(mode1) <- c("dec","inc")
      }
      if(any(data$diff>90)){
        mode2 <- as.data.frame(data$dec[data$diff>90])
        mode2$inc <- data$inc[data$diff>90]
        colnames(mode2) <- c("dec","inc")
      }
      if(exists("mode1")==TRUE) {fisher_M1 <- fisher(mode1)}
      if(exists("mode2")==TRUE) {fisher_M2 <- fisher(mode2)}
      if(plot==TRUE){
        if(exists("mode1")==TRUE){plot_a95(fisher_M1[1,1],fisher_M1[1,2],fisher_M1[1,3],
                                           on_plot = TRUE,symbol=symbol, col_d = col_d,
                                           col_u=col_u,col_l=col_l)}
        if(exists("mode2")==TRUE){plot_a95(fisher_M2[1,1],fisher_M2[1,2],fisher_M2[1,3],
                                           on_plot = TRUE,symbol=symbol, col_d = col_d,
                                           col_u=col_u,col_l=col_l)}
      }
      data_M12 <- common_DI(data)
      fisher_M12 <- fisher(data_M12)
      #plot text with results
      Dec <- round(fisher_M12[1,1],digits=2)
      Inc <- round(fisher_M12[1,2],digits=2)
      a <- round(fisher_M12[1,3],digits=2)
      N <- round(fisher_M12[1,4],digits=2)

      #creates table for Shiny
      S_results <- as.data.frame(matrix(ncol=6, nrow=3))
      colnames(S_results) <- c("dec", "inc", "a95", "N","R","k")


      if(any(data$diff<=90)) {
        S_results[1,] <- fisher_M1
      }
      if(any(data$diff>90)) {
        S_results[2,] <- fisher_M2
      }
      if(exists("fisher_M1")==TRUE | exists("fisher_M2")==TRUE) {
        S_results[3,] <- fisher_M12
      }
      rownames(S_results) <- c("Mode 1","Mode 2","All")
      S_results <- S_results[,-5]
      S_results <- na.omit(S_results)
      if(nrow(S_results)==2) S_results <- S_results[1,]

      return(S_results)
    }
    #modified ellips_plot function
    ellips_plot_S <- function(DI,lat=0,long=0, plot=TRUE, col_d="red",col_u="white",col_l="black",symbol="c"){
      #degrees to radians and vice versa
      d2r <- function(x) {x*(pi/180)}
      r2d <- function(x) {x*(180/pi)}
      #functions converting inc(x) and dec(y) into equal area
      a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
      a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
      data <- DI
      #cut lines with empty cells
      data <- na.omit(data)
      data <- data[,1:2]
      colnames(data) <- c("dec", "inc")
      #directions in Cartesian coordinates
      data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
      data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
      data$z <- sin(d2r(data$inc))
      #averaged Cartesian coordinates
      x_av <- mean(data$x)
      y_av <- mean(data$y)
      z_av <- mean(data$z)
      #elements of the distribution matrix
      T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                      sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                      sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
      #distribution matrix
      T <- matrix(T_elements,nrow=3, byrow=TRUE)
      #calculate and copy eigenvalues and vectors
      T_e <- eigen(T,symmetric = TRUE)
      T_vec <- T_e$vectors
      T_val <- T_e$values
      #calculate dec inc of max variance
      V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
      V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
      V1dec <- V1dec%%360
      #next  calculates difference between dec_inc and average
      data$Dec_aver <- rep(V1dec)
      data$Inc_aver <- rep(V1inc)
      data$delta <- abs(data$dec-data$Dec_aver)
      data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                              (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
      #Isolate modes
      if(any(data$diff<=90)){
        mode1 <- as.data.frame(data$dec[data$diff<=90])
        mode1$inc <- data$inc[data$diff<=90]
        colnames(mode1) <- c("dec","inc")
      }
      if(any(data$diff>90)){
        mode2 <- as.data.frame(data$dec[data$diff>90])
        mode2$inc <- data$inc[data$diff>90]
        colnames(mode2) <- c("dec","inc")
      }
      #calculate ellipses
      if(exists("mode1")==TRUE) {ellips_M1 <- ellips_DI(mode1, lat=lat, long=long)}
      if(exists("mode2")==TRUE) {ellips_M2 <- ellips_DI(mode2, lat=lat, long=long)}

      if(plot==TRUE){
        if(exists("mode1")==TRUE){generate_ellips(ellips_M1[1,1],ellips_M1[1,2],ellips_M1[1,3],ellips_M1[1,4],
                                                  on_plot = TRUE,symbol=symbol, col_d = col_d,
                                                  col_u=col_u,col_l=col_l)}
        if(exists("mode2")==TRUE){generate_ellips(ellips_M2[1,1],ellips_M2[1,2],ellips_M2[1,3],ellips_M2[1,4],
                                                  on_plot = TRUE,symbol=symbol, col_d = col_d,
                                                  col_u=col_u,col_l=col_l)}
      }
      data_M12 <- common_DI(data)
      ellips_M12 <- ellips_DI(data_M12,lat=lat, long=long)
      #plot text with results
      N <- ellips_M12[1,5]
      Dec <- round(ellips_M12[1,1],digits=2)
      Inc <- round(ellips_M12[1,2],digits=2)
      Delta_dec <- round(ellips_M12[1,3],digits=2)
      Delta_inc <- round(ellips_M12[1,4],digits=2)

      #set file for export in shiny
      S_result <- as.data.frame(matrix(ncol=5,nrow=3))
      colnames(S_result) <- c("dec", "inc","a95 dec","a95 inc","N")
      rownames(S_result) <- c("Mode 1","Mode 2", "All")

      if(any(data$diff<=90)) {
        S_result[1,] <- ellips_M1
      }
      if(any(data$diff>90)) {
        S_result[2,] <- ellips_M2
      }
      if(any(data$diff>90)) {
        S_result[3,] <- ellips_M12
      }
      S_result <- na.omit(S_result)
      if(nrow(S_result)==2) S_result <- S_result[1,]
      return(S_result)
    }

    output$directions <- renderPlot({
      #avoid errors if long and lat are missing
      req(input$lat)
      req(input$long)
      #create reactive file for stat
      F_stat <- reactiveValues(result=NULL)

      #import data
      DI <- fix_DI(input_file())

      #equal area function
      plot_dirs <- function(DI,Slat=input$lat,Slong=input$long,mode=input$mode,
                            colD=input$colD,colU=input$colU,sym=input$sym){
        #file with possible error message from inclination only routine
        inc_warn <- NULL
        if(mode==1){DI <- DI}
        if(mode==2){DI <- common_DI(DI)}
        if(mode==3){DI <- common_DI(DI,down = F)}
        #next just flip negative to positive or vice versa when required
        if(mode==4){
          for(i in 1:nrow(DI)){
            if(DI[i,2]<0){
              DI[i,1] <- (DI[i,1]+180)%%360
              DI[i,2] <- abs(DI[i,2])
            }
          }
        }else if(mode==5){
          for(i in 1:nrow(DI)){
            if(DI[i,2]>=0){
              DI[i,1] <- (DI[i,1]+180)%%360
              DI[i,2] <- -(DI[i,2])
            }
          }
        }
        #define colors Down-pointing
        if(colD==1) colD <- "black"
        if(colD==2) colD <- "blue"
        if(colD==3) colD <- "red"
        if(colD==4) colD <- "dark green"

        #define color Up-pointing
        if(colU==1) colU <- "white"
        if(colU==2) colU <- "cyan"
        if(colU==3) colU <- "pink"
        if(colU==4) colU <- "light green"

        #define symbol
        if(sym==1) sym <- "c"
        if(sym==2) sym <- "s"
        if(sym==3) sym <- "d"
        if(sym==4) sym <- "t"

        plot_DI(DI,col_d = colD,col_u = colU, symbol = sym)
        #plot statistic, with warning message if any from Arason+Levi2010 algorythm, or assign NULL in all other cases
        #inc warn must be created in the environment normally to cover the one created by PmagDiR
        if(input$fisher==2){
          assign("inc_warn",NULL, envir = .GlobalEnv)
          F_stat$result <- fisher_plot_S(DI)
        }else if(input$fisher==3){
          assign("inc_warn",NULL, envir = .GlobalEnv)
          F_stat$result <- ellips_plot_S(DI,lat = Slat,long = Slong)
        }else if(input$fisher==4){
          assign("inc_warn",NULL, envir = .GlobalEnv)
          F_stat$result <- PmagDiR::inc_plot(DI = DI,bimodal = F,print = F,export = F,save = F,Shiny = T)
        }else if(input$fisher==5){
          assign("inc_warn",NULL, envir = .GlobalEnv)
          F_stat$result <- PmagDiR::inc_plot(DI = DI,bimodal = T,print = F,export = F,save = F,Shiny = T)
        }else if(input$fisher==6){
          assign("inc_warn",NULL, envir = .GlobalEnv)
          F_stat$result <- PmagDiR::inc_plot(DI = DI,bimodal = F,print = F,export = F,save = F,arith_stat = T,Shiny = T)
        }else if(input$fisher==7){
          assign("inc_warn",NULL, envir = .GlobalEnv)
          F_stat$result <- PmagDiR::inc_plot(DI = DI,bimodal = T,print = F,export = F,save = F,arith_stat = T,Shiny = T)
        }else{
          assign("inc_warn",NULL, envir = .GlobalEnv)
        }
      }
      plot_dirs(DI)

      #record plot
      DirsPlot <- recordPlot()

      #create directions stats table
      output$stats <- renderTable({
        F_stat$result
      },rownames=T, digits=1)

      #write inc_only problem if any
      if(exists("inc_warn")==T){
        output$inc_warn <- renderText({
          inc_warn
        })
      }

      #export stat
      output$exportS <- downloadHandler(
        filename = function() {
          paste(input$fileN,"_", Sys.Date(), "_stat.csv", sep="")
        },
        content = function(file) {
          write.csv(round(F_stat$result, digits = 2), file)
        }
      )

      #export equal_area graph
      output$exportG <- downloadHandler(
        filename = function() {
          paste(input$fileN,"_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 9,height = 9)
          replayPlot(DirsPlot)
          dev.off()
        }
      )

      #export DI
      output$exportDI <- downloadHandler(
        filename = function() {
          paste(input$fileN,"_", Sys.Date(), "_directions.csv", sep="")
        },
        content = function(file) {
          DI <- fix_DI(input_file())
          if(input$mode==1){DI <- DI}
          if(input$mode==2){DI <- common_DI(DI)}
          if(input$mode==3){DI <- common_DI(DI,down = F)}
          DI <- na.omit(DI)
          write.csv(round(DI, digits=2),row.names = F, file)
        }
      )
    },width = 700,height = 700)
    ############ END OF EQUAL AREA MODULE

    ############ REVERSAL TEST MODULE
    #create function that reversal test and save statistics and graph
    revtest_plot <- eventReactive(input$revgo,{

      #main revtest function adapted for Shiny
      revtest_S <- function(DI,nb){
        #fucnctions deg to rads and vice versa
        d2r <- function(x) {x*(pi/180)}
        r2d <- function(x) {x*(180/pi)}

        data <- DI[,1:2]
        data <- na.omit(data)
        colnames(data) <- c("dec", "inc")

        #directions in Cartesian coordinates
        data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
        data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
        data$z <- sin(d2r(data$inc))

        #averaged Cartesian coordinates
        x_av <- mean(data$x)
        y_av <- mean(data$y)
        z_av <- mean(data$z)

        #elements of the distribution matrix
        T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                        sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                        sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))

        #distribution matrix
        T <- matrix(T_elements,nrow=3, byrow=TRUE)

        #calculate and copy eigenvalues and vectors
        T_e <- eigen(T,symmetric = TRUE)
        T_vec <- T_e$vectors
        T_val <- T_e$values

        #calculate dec inc of max variance
        V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
        V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
        V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)

        #flip V1 if negative
        V1dec <- ifelse(V1inc<0,ifelse((V1dec+180)>360,V1dec-180,V1dec+180),V1dec)
        V1inc <- ifelse(V1inc<0,-V1inc,V1inc)


        #next  calculates difference between dec_inc and average
        data$Dec_aver <- rep(V1dec)
        data$Inc_aver <- rep(V1inc)
        data$delta <- abs(data$dec-data$Dec_aver)
        data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                                (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
        #Isolate modes
        m1ind <- as.numeric(which(data$diff<=90), arr.ind = TRUE)
        m2ind <- as.numeric(which(data$diff>90), arr.ind = TRUE)

        #terminate if distribution is not bimodal
        if(length(m2ind)<1) stop("
DISTRIBUTION NOT BIMODAL")
        mode1 <- data[m1ind,1:2]
        mode2 <- data[m2ind,1:2]

        #flip mode 2 same as mode 1
        mode2 <- flip_DI(mode2)
        nb <- nb
        n1 <- 0
        n2 <- 0

        mode1B <- as.data.frame(matrix(ncol=3, nrow=0))
        mode2B <- as.data.frame(matrix(ncol=3, nrow=0))

        #simulate pseudosamples of mode 1
        repeat{
          n1 <- n1+1
          mode1B_p <- as.data.frame(matrix(ncol=3, nrow=1))
          Bdata <- boots_DI(mode1)
          Bdata$x <- cos(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
          Bdata$y <- sin(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
          Bdata$z <- sin(d2r(Bdata$inc))
          mode1B_p[1,1] <- mean(Bdata$x)
          mode1B_p[1,2] <- mean(Bdata$y)
          mode1B_p[1,3] <- mean(Bdata$z)
          mode1B <- rbind(mode1B,mode1B_p)
          #function that update the progress bar of shiny
          updateProgressBar(
            id="mode1",
            value=n1,total=nb,
          )
          if(n1==nb) break
        }
        colnames(mode1B) <- c("x","y","z")
        mode1B$dec <- r2d(atan2(mode1B$y,mode1B$x))
        mode1B$dec <- mode1B$dec%%360
        mode1B$inc <- r2d(asin(mode1B$z))
        #simulate pseudosamples of mode 2
        repeat{
          n2 <- n2+1
          mode2B_p <- as.data.frame(matrix(ncol=3, nrow=1))
          Bdata <- boots_DI(mode2)
          Bdata$x <- cos(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
          Bdata$y <- sin(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
          Bdata$z <- sin(d2r(Bdata$inc))
          mode2B_p[1,1] <- mean(Bdata$x)
          mode2B_p[1,2] <- mean(Bdata$y)
          mode2B_p[1,3] <- mean(Bdata$z)
          mode2B <- rbind(mode2B,mode2B_p)
          updateProgressBar(
            id="mode2",
            value=n2,total=nb,
          )
          if(n2==nb) break
        }
        colnames(mode2B) <- c("x","y","z")
        mode2B$dec <- r2d(atan2(mode2B$y,mode2B$x))
        mode2B$dec <- ifelse(mode2B$dec<0,mode2B$dec+360,mode2B$dec)
        mode2B$inc <- r2d(asin(mode2B$z))

        #isolate components of models
        B1x <- sort(mode1B[,1])
        B1y <- sort(mode1B[,2])
        B1z <- sort(mode1B[,3])
        B2x <- sort(mode2B[,1])
        B2y <- sort(mode2B[,2])
        B2z <- sort(mode2B[,3])

        #define low and high boostrapped margins
        confn <- 0.95
        num <- round((nb*(1-confn))/2,digits=0)
        Lconf <- num
        Uconf <- nb-num
        B1x_l <- c(B1x[Lconf],B1x[Uconf])
        B2x_l <- c(B2x[Lconf],B2x[Uconf])
        B1y_l <- c(B1y[Lconf],B1y[Uconf])
        B2y_l <- c(B2y[Lconf],B2y[Uconf])
        B1z_l <- c(B1z[Lconf],B1z[Uconf])
        B2z_l <- c(B2z[Lconf],B2z[Uconf])

        #max and min values for graphs
        xmax <- round(max(c(B1x,B2x)), digits=1)+0.05
        xmin <- round(min(c(B1x,B2x)),digits=1)-0.05
        ymax <- round(max(c(B1y,B2y)), digits=1)+0.05
        ymin <- round(min(c(B1y,B2y)),digits=1)-0.05
        zmax <- round(max(c(B1z,B2z)), digits=1)+0.05
        zmin <- round(min(c(B1z,B2z)),digits=1)-0.05

        #function that extract intervals and counts from hist function and make cumulative curve
        cumulative_curve <- function(x){
          h <- hist(x, breaks=50,plot = FALSE)
          cnts <- h[["counts"]]
          t <- length(cnts)
          new_c <- as.data.frame(matrix(ncol=1,nrow = 1))
          for(i in 1:t){
            if(i==1){new_cp <- cnts[1]}
            if(i>1) {new_cp <- new_c[i,1]+cnts[i-1]}
            new_c <- rbind(new_c,new_cp)
          }
          new_c <- na.omit(new_c)
          breaks <- as.data.frame(h[["mids"]])
          cumul <- cbind(breaks,new_c)
          colnames(cumul) <- c("breaks","counts")
          cumul$counts <- cumul$counts/nb
          return(cumul)
        }
        cu1x <- cumulative_curve(B1x)
        cu2x <- cumulative_curve(B2x)
        cu1y <- cumulative_curve(B1y)
        cu2y <- cumulative_curve(B2y)
        cu1z <- cumulative_curve(B1z)
        cu2z <- cumulative_curve(B2z)
        text1 <- "Pseudosample means"
        text2 <- "Normalized cumulative distributions"

        #clean screen to avoid figure over figure
        par(fig=c(0,1,0,1))
        plot(0, xlim=c(0,1), ylim=c(0,1),
             xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

        #plot title for equal area
        par(fig=c(0,0.55,0.4,1))
        plot(NA, xlim=c(0,1), ylim=c(0,1),
             xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
        text(x=0.5,y=0.8,substitute(paste(bold(text1))),cex=1.5,pos=3)

        #plot title for cumulative distributions
        par(fig=c(0.55,1,0.4,1),new=TRUE)
        plot(NA, xlim=c(0,1), ylim=c(0,1),
             xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
        text(x=0.5,y=0.8,substitute(paste(bold(text2))),cex=1.5, pos=3)

        #plot equal area
        #pseudosamples mean
        par(fig=c(0,0.55,0,0.9), new=TRUE)
        plot_DI(mode1B[,4:5],col_d = rgb(1,0,0,0.20),col_u=rgb(1,0.75,1,0.30),col_ext = NA)
        plot_DI(mode2B[,4:5],col_d =rgb(0,0,1,0.20),col_u=rgb(0,1,1,0.30),col_ext = NA,on_plot = TRUE)

        #plot cumulative distributions
        par(fig=c(0.55,1,0.52,0.82),new=TRUE)
        #plot x
        plot(0,type="n",xlim=c(xmin,xmax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
        title(xlab = "x axis", line=2, cex=1.2)
        rect(xleft = B1x_l[1],xright=B1x_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
        rect(xleft = B2x_l[1],xright=B2x_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
        points(cu1x, type="l", lwd=2, col="red")
        points(cu2x, type="l", lwd=2, col="blue")

        par(fig=c(0.55,1,0.31,0.61),new=TRUE)
        #plot y
        plot(0,type="n",xlim=c(ymin,ymax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
        title(xlab = "y axis", line=2, cex=1.2)
        rect(xleft = B1y_l[1],xright=B1y_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
        rect(xleft = B2y_l[1],xright=B2y_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
        points(cu1y, type="l", lwd=2, col="red")
        points(cu2y, type="l", lwd=2, col="blue")

        par(fig=c(0.55,1,0.1,0.4),new=TRUE)
        #plot z
        plot(0,type="n",xlim=c(zmin,zmax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
        title(xlab = "z axis", line=2, cex=1.2)
        rect(xleft = B1z_l[1],xright=B1z_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
        rect(xleft = B2z_l[1],xright=B2z_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
        points(cu1z, type="l", lwd=2, col="red")
        points(cu2z, type="l", lwd=2, col="blue")
        #reset screen
        par(fig=c(0,1,0,1))

        revtest_results <- list(0)
        #export
        boot_stat <- as.data.frame(matrix(c(B1x_l[1],B1x_l[2],B2x_l[1],B2x_l[2],
                                            B1y_l[1],B1y_l[2],B2y_l[1],B2y_l[2],
                                            B1z_l[1],B1z_l[2],B2z_l[1],B2z_l[2]),nrow = 3,byrow = TRUE))
        rownames(boot_stat) <- c("x","y","z")
        colnames(boot_stat) <- c("M1_L","M1_H","M2_L","M2_H")
        revtest_results[[1]] <- boot_stat
        return(revtest_results)
      }

      #reversal test
      revtest_funct <- function(){
        DI <- fix_DI(input_file())
        revtest_S(DI,nb=input$revnb)
      }

      #save statistic and makes plot
      revdat <- revtest_funct()
      #record plot
      revPlot <- recordPlot()
      #prepare table for display
      output$revstat <- renderTable({revdat[[1]]},rownames = T, digits = 2,
                                    caption="M1,2 = Mode 1 and 2. L and H = lower and higher confidence limit")
      #Export reversal test graphic
      output$revexpG <- downloadHandler(
        filename = function() {
          paste(input$fileN_RT,"_revtest_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 15,height = 10)
          replayPlot(revPlot)
          dev.off()
        }
      )
      #Export reversal test stats
      output$revexpS <- downloadHandler(
        filename = function() {
          paste(input$fileN_RT,"_revtest_", Sys.Date(), "_stat.csv", sep="")
        },
        content = function(file) {
          write.csv(round(revdat[[1]],digits = 3), file)
        }
      )


    })
    #execute reversal test
    output$revtest <- renderPlot({
      revtest_plot()
    },width = 1000,height = 700)
    ############ END OF REVERSAL TEST MODULE

    ############ DISTRIBUTION SHAPE MODULE
    #reactive result file
    EI_boot <- reactiveValues(result=NULL)

    #main  function adapted for Shiny
    EI_boot_S <- function(DI,boot=input$EIyesnoboot,nb=1000,conf=95) {
      data <- DI
      data <- na.omit(data)
      data <- data[,1:2]
      colnames(data) <- c("dec", "inc")
      Inc_E_real <- PmagDiR::inc_E_finder(data)
      Inc_E_real$V1inc <- abs(Inc_E_real$V1inc)
      Inc_E <- as.data.frame(matrix(ncol=3,nrow=0))

      #plot frame
      par(fig=c(0,1,0,1), new= FALSE)
      y_up <- ifelse(Inc_E_real$E>3.5, ceiling(Inc_E_real$E),3.5)
      plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
           xlab="Inclination (°)", ylab="Elongation", cex.lab=1.5)

      #Plot real data E-I couple and values
      points(x=Inc_E_real$V1inc,y=Inc_E_real$E,pch=21,
             col="black", bg="red", cex=1.5)

      #perform bootstrap only if requested
      if(boot==2){
        for (i in 1:nb){
          dataprov <- boots_DI(data)
          I_E_Ed <- inc_E_finder(dataprov)
          I_E_Ed$V1inc <- abs(I_E_Ed$V1inc)
          Inc_E <- rbind(Inc_E, I_E_Ed)
          #function that update the progress bar of shiny
          updateProgressBar(
            id="EIbootstrap",
            value=i,total=nb,
          )
          if(i==nb) break
        }
        colnames(Inc_E) <- c("Inc","E","E_dec")
        Inc_E <- Inc_E[order(Inc_E$E),]
        Inc_E_bk <- Inc_E
        Inc_E <- Inc_E_bk

        confn <- conf/100
        num <- round((nb*(1-confn))/2,digits=0)
        Lconf <- num
        Uconf <- nb-num
        Inc_E <- Inc_E[Lconf:Uconf,]

        #Plot Bootstrapped data
        points(x=Inc_E$Inc,
               y=Inc_E$E,
               pch=16,
               col=rgb(0, 0, 1, 0.05),
               cex=0.8)

        #Re-plot real data E-I couple and values
        points(x=Inc_E_real$V1inc,y=Inc_E_real$E,pch=21,
               col="black", bg="white", cex=1.5)
      }

      #plot tk03.GAD model E-I
      x <- 0:90
      y <- tk03(x)
      points(x=x, y= y, type= "l", col="blue", lwd=3)

      N <- nrow(data)
      Inc <- format(round(Inc_E_real$V1inc,1),nsmall=1)
      Ecut <- format(round(Inc_E_real$E,2),nsmall=2)
      V2 <- format(round(Inc_E_real$DV1V2,1),nsmall=1)

      text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
      text(x=0, y=3.2,pos=4,text, cex= 1.2)

      if(boot==2){
        #plot histogram of Dec top right
        par(fig=c(0.6,1,0.56,0.99), new=TRUE)

        hist(Inc_E$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
             breaks = 90, xlab = "", ylab = "",
             main="", cex.axis=0.8, col= "blue", border="blue")
        #plot labels closer than standard to axes
        title(xlab = "Edec(°)", line=1.9, cex=0.5)
        title(ylab = "Frequency", line=2, cex=0.3)

        #plot confidence margin of declination if between -45 and 45
        if(Inc_E_real$DV1V2>-50 && Inc_E_real$DV1V2<50){
          Inc_E_Edec <- Inc_E_bk
          Inc_E_Edec <- Inc_E_Edec[order(Inc_E_Edec[,3]),]
          Inc_E_Edec <- Inc_E_Edec[Lconf:Uconf,]
          low_dec <- Inc_E_Edec[1,3]
          up_dec <- Inc_E_Edec[length(Inc_E_Edec[,3]),3]
          abline(v=low_dec,lwd=1,lty=2)
          abline(v=up_dec,lwd=1,lty=2)
        }
        #calculate low and high inclination error
        Inc_E <- Inc_E_bk
        Inc_E <- Inc_E[order(Inc_E[,1]),]
        Inc_E <- Inc_E[Lconf:Uconf,]
        low_inc <- Inc_E[1,1]
        up_inc <- Inc_E[length(Inc_E[,1]),1]
        EI_boot_tab <- as.data.frame(matrix(ncol= 3, nrow=3))
        colnames(EI_boot_tab) <- c("Mean","Low","High")
        rownames(EI_boot_tab) <- c("Inc","Elong","E_dec")
        EI_boot_tab[1,1] <- Inc
        EI_boot_tab[1,2] <- round(low_inc,digits = 1)
        EI_boot_tab[1,3] <- round(up_inc,digits = 1)
        EI_boot_tab[2,1] <- Ecut
        EI_boot_tab[2,2] <- round(min(Inc_E$E),digits= 2)
        EI_boot_tab[2,3] <- round(max(Inc_E$E),digits= 2)
        EI_boot_tab[3,1] <- V2
        if(Inc_E_real$DV1V2>-50 && Inc_E_real$DV1V2<50){
          EI_boot_tab[3,2] <- round(low_dec,digits= 2)
          EI_boot_tab[3,3] <- round(up_dec,digits=2)
        }else{
          EI_boot_tab[3,2] <- ""
          EI_boot_tab[3,3] <- ""
        }
        #reset screen
        return(EI_boot_tab)
      }
    }

    #reactive function
    EI_boot_plot <- eventReactive(input$EIbootgo,{

      #perform routine
      EI_boot_funct <- function(){
        DI <- fix_DI(input_file(),coord=2)
        EI_boot_S(DI,nb=input$EIbootnb,conf = input$EIconf)
      }

      #save statistic and makes plot
      EI_boot$result <- EI_boot_funct()
      #prepare table for display
      output$EIbootStat <- renderTable({EI_boot$result},rownames = T, digits = 2,
                                       caption="Inc= inclination, Elong= elongation, E_dec= declination of elongation")
    })

    #execute EIboot test
    output$EIboot <- renderPlot({
      EI_boot_plot()
      #record plot
      EI_bootPlot <- recordPlot()
      #Export graphic
      output$EIbootG <- downloadHandler(
        filename = function() {
          paste(input$fileN_EI,"_EIboot_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 12,height = 11)
          replayPlot(EI_bootPlot)
          dev.off()
        }
      )
      #Export stats
      output$EIbootS <- downloadHandler(
        filename = function() {
          paste(input$fileN_EI,"_EIboot_", Sys.Date(), "_stat.csv", sep="")
        },
        content = function(file) {
          write.csv(EI_boot$result, file)
        }
      )
    }, width = 800,height = 700)
    ############ END OF DISTRIBUTION SHAPE MODULE

    ############ INCLINATION FLATTENING MODULE
    #Function that correct inclination shallowing after tk03.GAD model
    #ffind function is external cause it is used by the VGPs module
    ffind_boot_S <- function(DI,confidence=95,nb=1000, bootstrap= 1, f_increment=0.01,strat_pos=F) {
      data <- DI[,1:2]
      data <- na.omit(data)
      N <- length(data[,1])
      colnames(data) <- c("dec", "inc")
      #calculate E-I of real data
      Inc_E_R <- inc_E_finder(data)
      Inc_E_R$V1inc <- abs(Inc_E_R$V1inc)
      Inc <- round(Inc_E_R$V1inc, digits=1)
      Ecut <- round(Inc_E_R$E, digits=2)
      Edec <- round(Inc_E_R$DV1V2, digits=1)

      #calculate E-I correction sequence of real data
      Seq_I_E_R <- ffind(data,f_inc = 0.0005)
      colnames(Seq_I_E_R) <- c("V1inc","E","DV1V2","f")
      alert <- ifelse(length(Seq_I_E_R$V1inc)==1,"y","n")
      Ffinal <- round(Seq_I_E_R[length(Seq_I_E_R$f),4], digits=2)
      Inc_f <- round(Seq_I_E_R[length(Seq_I_E_R$V1inc),1], digits=1)
      Efinal <- round(Seq_I_E_R[length(Seq_I_E_R$E),2], digits=2)
      Edec_f <- round(Seq_I_E_R[length(Seq_I_E_R$DV1V2),3], digits=1)
      f <- min(Seq_I_E_R$f)
      if(alert=="y") f <- 1
      unf_data <- PmagDiR::unflat_DI(data,f)
      colnames(unf_data) <- c("dec","inc")

      #plot frame
      par(fig=c(0,1,0,1), new= FALSE)
      plot(NA, xlim= c(0,90), ylim= c(1,3.5), xaxp= c(0,90, 9),
           xlab="Inclination (°)", ylab="Elongation",cex.lab=1.5)

      #plot tk03.GAD model E-I
      x <- 0:90
      y <- tk03(x)
      points(x=x, y= y, type= "l", col="blue", lwd=3)

      #Plot real data E-I and correction of real data
      points(x=Seq_I_E_R$V1inc, y= Seq_I_E_R$E, type= "l", col="red", lwd=3)
      points(x=Inc_E_R$V1inc,y=Inc_E_R$E,pch=21,
             col="black", bg="blue", cex=1.5)
      points(x=Inc_f,y=Efinal,
             pch=21, col="black", bg="red", cex=1.5)

      text <- paste("N:", N, "
Inc:", Inc,"
E:", Ecut,"
Edec:",Edec)

      text2 <- paste("f:", Ffinal, "
Inc_Unfl:", Inc_f, "
E_Unfl:", Efinal, "
Edec_Unfl:", Edec_f)

      text(x=0, y=3.2,pos=4,text, cex= 1.2)
      text(x=20, y=3.2, pos=4, text2, cex=1.2)
      if (alert=="y"){
        text3 <- "Distribution not flattened"
        text(x=0, y=3, pos=4,text3,cex=1)
      }
      #create files for initial and final readings for the histograms
      init_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
      final_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
      colnames(init_E_I) <- c("Inc","E","E_dec")
      colnames(final_E_I)<- c("Inc","E","E_dec")
      if(alert=="y") par(fig=c(0,1,0,1))
      if(alert=="y") stop("
DISTRIBUTION NOT FLATTENED.")
      #creates file for fisher of bootstrapped pseudosamples
      b_fisher <- data.frame(matrix(ncol=6, nrow=0))
      n <- 0
      par(fig=c(0,1,0,1), new=TRUE)
      plot(NA, xlim= c(0,90), ylim= c(1,3.5), xaxt="n",yaxt="n",
           xlab="", ylab="", axes=FALSE)
      if(bootstrap==2){
        repeat {
          n <- n+1
          Seq_I_E_B <- as.data.frame(matrix(ncol=3,nrow=0))
          dataprov <- boots_DI(data)
          Seq_I_E_B <- ffind(dataprov, f_inc = f_increment)
          #calculate fisher of bootstrapped flattened dataset and paste on file
          b_fisher_temp <- fisher(unflat_DI(common_DI(dataprov),f = min(Seq_I_E_B[,4])))
          b_fisher <- rbind(b_fisher,b_fisher_temp)
          #plot bootstrapped lines
          points(x=Seq_I_E_B$V1inc, y= Seq_I_E_B$E,
                 type= "l", col=rgb(1, 0, 0, 0.08), lwd=1.2)
          i_E_I <- Seq_I_E_B[1,]
          f_E_I <- Seq_I_E_B[length(Seq_I_E_B[,1]),]
          colnames(i_E_I) <- c("Inc","E","E_dec")
          colnames(f_E_I)<- c("Inc","E","E_dec")

          #isolate initial and final readings for histograms
          init_E_I <- rbind(init_E_I,i_E_I)
          final_E_I <- rbind(final_E_I,f_E_I)
          init_E_I <- na.omit(init_E_I)
          final_E_I <- na.omit(final_E_I)
          finaln <- length(final_E_I[,1])
          updateProgressBar(
            id="ffindbootstrap",
            value=finaln,total=nb,
          )

          if(finaln==nb) {
            output$validboots <- renderText({paste("Total number of simulations:",n)})
            break
          }
        }

        #replot real data with different color
        points(x=Seq_I_E_R$V1inc, y= Seq_I_E_R$E, type= "l", col="yellow", lwd=3)
        points(x=Inc_E_R$V1inc,y=Inc_E_R$E,pch=21,
               col="black", bg="blue", cex=1.5)
        points(x=x, y= y, type= "l", col="blue", lwd=3)
        points(x=Inc_f,y=Efinal,
               pch=21, col="black", bg="red", cex=1.5)


        #replot results in case covered by boostrapps
        text(x=0, y=3.2,pos=4,text, cex= 1.2)
        text(x=20, y=3.2, pos=4, text2, cex=1.2)

        colnames(init_E_I) <- c("Inc","E","E_dec","f")
        colnames(final_E_I)<- c("Inc","E","E_dec","f")
        final_E_I <- final_E_I[order(final_E_I$Inc),] #order final results by inclination
        final_E_Ibk <- final_E_I

        conf <- confidence/100
        num <- round((nb*(1-conf))/2,digits=0)
        Lconf <- num
        Uconf <- nb-num
        final_E_I <- final_E_I[Lconf:Uconf,]    #cut bootstrapped results for 95% confidence

        #draw two lines for 95% confidence margin
        arrows(x0=final_E_I[1,1],x1=final_E_I[1,1],
               y0=1,y1=final_E_I[1,2], length = 0,lty=2)

        arrows(x0=final_E_I[length(final_E_I$Inc),1],
               x1=final_E_I[length(final_E_I$Inc),1],
               y0= 1, y1= final_E_I[length(final_E_I$Inc),2],
               length = 0,lty=2)
        Inc_l95 <- round(final_E_I[1,1], digits= 1)
        Inc_u95 <- round(final_E_I[length(final_E_I$Inc),1], digits=1)

        text(x=final_E_I[1,1], y=1, pos=2, Inc_l95,cex=1.2)
        text(x=final_E_I[length(final_E_I$Inc),1], y=1, pos=4, Inc_u95,cex=1.2)


        #plot histogram of E_declination with respect V1 before and after correction
        par(fig=c(0.6,1,0.56,0.99), new=TRUE)
        hist(init_E_I$E_dec, xlim=c(-90,90), breaks= 90,
             axes=FALSE,xlab="",ylab="",col="blue", border="blue", main="")
        #plot lables closer than standard to axes
        title(xlab = "Edec(°)", line=1.9, cex=0.2)
        title(ylab = "Frequency", line=1.9,cex=0.2)
        #after
        par(fig=c(0.6,1,0.56,0.99), new=TRUE)
        hist(final_E_Ibk$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
             breaks = 90, xlab = "", ylab = "",
             main="", cex.axis=0.8,col="red",border ="red")

        #recalculate boostrtapped confidence for Edec for plotting confidence margin
        final_E_I_Edec <- final_E_Ibk
        colnames(final_E_I_Edec) <- c("Inc","E","E_dec")
        final_E_I_Edec <- final_E_I_Edec[order(final_E_I_Edec$E_dec),]
        final_E_I_Edec <- final_E_I_Edec[Lconf:Uconf,]
        abline(v=final_E_I_Edec[1,3],lwd=1, lty=2)
        abline(v=final_E_I_Edec[length(final_E_I_Edec[,3]),3],lwd=1, lty=2)

        #plot histogram of inclination before and after correction
        par(fig=c(0.6,1,0.25,0.68), new=TRUE)
        hist(init_E_I$Inc,xlim=c(0,90), breaks=90,
             axes=FALSE,xlab="",ylab="",col="blue",border="blue", main="")
        par(fig=c(0.6,1,0.25,0.68), new=TRUE)
        hist(final_E_Ibk$Inc, xlim=c(0,90), xaxp=c(0,90,6),
             breaks=90,xlab = "", ylab = "",
             main="",cex.axis=0.8,col="red",border ="red")
        abline(v=final_E_I[1,1],lwd=1,lty=2)
        abline(v=final_E_I[length(final_E_I[,1]),1],lwd=1,lty=2)
        title(xlab = "Inc(°)", line=1.9, cex=0.2)
        title(ylab = "Frequency", line=1.9, cex=0.2)
      }

      ffind_boot_stat <- data.frame(matrix(ncol=3,nrow = 5))
      colnames(ffind_boot_stat) <- c("Mean","Low","High")
      rownames(ffind_boot_stat) <- c("Inc","E_dec","f","Elong","N")
      ffind_boot_stat[1,1] <- Inc_f
      ffind_boot_stat[1,2] <- ifelse(bootstrap==2,Inc_l95,"")
      ffind_boot_stat[1,3] <- ifelse(bootstrap==2,Inc_u95,"")
      ffind_boot_stat[2,1] <- Edec_f
      ffind_boot_stat[2,2] <- ifelse(bootstrap==2,round(final_E_I_Edec[1,3],digits = 2),"")
      ffind_boot_stat[2,3] <- ifelse(bootstrap==2,round(final_E_I_Edec[nrow(final_E_I_Edec),3],digits=2),"")
      ffind_boot_stat[3,1] <- round(f, digits=2)
      ffind_boot_stat[3,2:3] <- t(c("",""))
      ffind_boot_stat[4,1] <- Efinal
      ffind_boot_stat[4,2:3] <- t(c("",""))
      ffind_boot_stat[5,1] <- N
      ffind_boot_stat[5,2:3] <- t(c("",""))
      #create result list for export
      ffind_boot_result <- list(0)
      ffind_boot_result[[1]] <- unf_data
      ffind_boot_result[[2]] <- ffind_boot_stat
      ffind_boot_result[[3]] <- b_fisher
      return(ffind_boot_result)

      par(fig=c(0,1,0,1))
    }
    #create reactive file
    ffind <- reactiveValues(results=NULL)
    ffind_boot_plot <- eventReactive(input$ffindgo,{
      #main function adapted
      #perform routine
      ffind_boot_funct <- function(){
        DI <- fix_DI(input_file(),coord=2)
        ffind_boot_S(DI,nb=input$ffindboot,bootstrap=input$ffindyesnoboot)
      }
      #save statistic and makes plot
      ffind$results <- ffind_boot_funct()
      #prepare and plot table
      output$ffindStat <- renderTable({ffind$results[[2]][1:2,]},rownames = T,
                                      caption="Inc= inclination, E_dec= declination of elongation. Exported file includes f, N, Elongation.")
      #plot corrected directions, always tilt corrected coords
      output$directions_FF <- renderPlot({
        #equal area function
        plot_EA <- function(dirs){
          req(input$lat)
          req(input$long)
          Slat <- input$lat
          Slong <- input$long
          DI <- dirs
          if(input$mode_FF==1){DI <- DI}
          if(input$mode_FF==2){DI <- common_DI(DI)}
          if(input$mode_FF==3){DI <- common_DI(DI,down = F)}
          #next just flip negative to positive or vice versa when required
          if(input$mode_FF==4){
            for(i in 1:nrow(DI)){
              if(DI[i,2]<0){
                DI[i,1] <- (DI[i,1]+180)%%360
                DI[i,2] <- abs(DI[i,2])
              }
            }
          }else if(input$mode_FF==5){
            for(i in 1:nrow(DI)){
              if(DI[i,2]>=0){
                DI[i,1] <- (DI[i,1]+180)%%360
                DI[i,2] <- -(DI[i,2])
              }
            }
          }
          #define colors Down-pointing
          if(input$colD_FF==1) colD <- "black"
          if(input$colD_FF==2) colD <- "blue"
          if(input$colD_FF==3) colD <- "red"
          if(input$colD_FF==4) colD <- "dark green"

          #define color Up-pointing
          if(input$colU_FF==1) colU <- "white"
          if(input$colU_FF==2) colU <- "cyan"
          if(input$colU_FF==3) colU <- "pink"
          if(input$colU_FF==4) colU <- "light green"

          #define symbol
          if(input$sym_FF==1) sym <- "c"
          if(input$sym_FF==2) sym <- "s"
          if(input$sym_FF==3) sym <- "d"
          if(input$sym_FF==4) sym <- "t"

          plot_DI(DI,col_d = colD,col_u = colU, symbol = sym)
          if(input$fisher_FF==2){
            fisher_plot_S(DI)
          }else if(input$fisher_FF==3){
            ellips_plot_S(DI,lat = Slat,long = Slong)
          }
        }
        plot_EA(ffind$results[[1]])
      },width = 700,height = 700)
      #Export graphic  equal area
      output$ffindG_FF <- downloadHandler(
        filename = function() {
          paste(input$fileN_FF,"_unflatdirs_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 9,height = 9)
          plot_EA <- function(dirs){
            req(input$lat)
            req(input$long)
            Slat <- input$lat
            Slong <- input$long
            DI <- dirs
            if(input$mode_FF==1){DI <- DI}
            if(input$mode_FF==2){DI <- common_DI(DI)}
            if(input$mode_FF==3){DI <- common_DI(DI,down = F)}
            #define colors Down-pointing
            if(input$colD_FF==1) colD <- "black"
            if(input$colD_FF==2) colD <- "blue"
            if(input$colD_FF==3) colD <- "red"
            if(input$colD_FF==4) colD <- "dark green"

            #define color Up-pointing
            if(input$colU_FF==1) colU <- "white"
            if(input$colU_FF==2) colU <- "cyan"
            if(input$colU_FF==3) colU <- "pink"
            if(input$colU_FF==4) colU <- "light green"

            #define symbol
            if(input$sym_FF==1) sym <- "c"
            if(input$sym_FF==2) sym <- "s"
            if(input$sym_FF==3) sym <- "d"
            if(input$sym_FF==4) sym <- "t"

            plot_DI(DI,col_d = colD,col_u = colU, symbol = sym)
            if(input$fisher_FF==2){
              fisher_plot_S(DI)
            }else if(input$fisher_FF==3){
              ellips_plot_S(DI,lat = Slat,long = Slong)
            }
          }
          plot_EA(ffind$results[[1]])
          dev.off()
        }
      )
      #export unflattened directions
      output$ffindDI_FF <- downloadHandler(
        filename = function() {
          paste(input$fileN_FF,"_unflatdirs_", Sys.Date(),".csv", sep="")
        },
        content = function(file) {
          if(input$mode_FF==1){DI_FF <- ffind$results[[1]]}
          else if(input$mode_FF==2){DI_FF <- common_DI(ffind$results[[1]])}
          else if(input$mode_FF==3){DI_FF <- common_DI(ffind$results[[1]],down = F)}
          write.csv(round(DI_FF, digits = 2),file = file, row.names = F)
        })
      #expor unfl direction statistics
      output$ffindS_FF <- downloadHandler(
        filename = function(){
          paste(input$fileN_FF,"_unflatDIstat_", Sys.Date(),".csv", sep="")
        },
        content = function(file){
          req(input$lat)
          req(input$long)
          Slat <- input$lat
          Slong <- input$long
          DI <- ffind$results[[1]]
          if(input$fisher_FF==2){
            F_stat_FF <- fisher_plot_S(DI,plot = F)
          }else if(input$fisher_FF==3){
            F_stat_FF <- ellips_plot_S(DI,lat = Slat,long = Slong,plot=F)
          }else{F_stat_FF <- NULL}
          write.csv(round(F_stat_FF, digits=2),file = file)})

      #print statistic table
      output$stats_FF <- renderTable({
        req(input$lat)
        req(input$long)
        Slat <- input$lat
        Slong <- input$long
        DI <- ffind$results[[1]]
        if(input$fisher_FF==2){
          F_stat_FF <- fisher_plot_S(DI,plot = F)
        }else if(input$fisher_FF==3){
          F_stat_FF <- ellips_plot_S(DI,lat = Slat,long = Slong,plot=F)
        }else{F_stat_FF <- NULL}
        F_stat_FF
      },rownames=T, digits=1)
    })
    #execute ffind_boot test
    output$ffindgraph <- renderPlot({
      ffind_boot_plot()
      ffindPlot <- recordPlot()
      #Export graphic
      output$ffindG <- downloadHandler(
        filename = function() {
          paste(input$fileN_FF,"_ffind_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 12,height = 11)
          replayPlot(ffindPlot)
          dev.off()
        }
      )
      output$ffindS <- downloadHandler(
        filename = function() {
          paste(input$fileN_FF,"_ffind_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(ffind$results[[2]],file)
        }
      )
    },width = 800,height = 700)
    ############ END OF INCLINATION FLATTENING MODULE

    ############ VIRTUAL GEOMAGNETIC POLES MODULE
    #modified VGP plot function ; VGPint define if VGP calculation are coming from VGP screen 1,2, or 3
    plot_VGP_S <- function(VGP,lat=90,long=0,grid=30,plot_vgp=TRUE, symbol="c",cex=1,
                           col_sym_out="black", col="black",col_f="cyan",col_boot=rgb(1,0,0,0.15),
                           on_plot=FALSE,auto_cent=TRUE,coast=FALSE, title="",save=TRUE,A95=FALSE,B95=FALSE,nb=1000,VGPint=1,plot=T){
      #functions converting degree and radians
      d2r <- function(x) {x*(pi/180)}
      r2d <- function(x) {x*(180/pi)}

      #functions converting long & lat to xy
      c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
      c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
      #cut is cosin of c, when negative is behind projections, needs to be cut
      cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
      #manipulate data
      colnames(VGP) <- c("lon","lat")
      vgpsN <- common_DI(VGP,down = ifelse(mean(VGP$lat)<0,FALSE,TRUE))
      PPole <- fisher(vgpsN)
      #fix point of view
      if(auto_cent==FALSE){
        #center of proj is Lon0 & Lat0
        lon0 <- long
        lat0 <- lat
      }else{
        lon0 <- PPole[1,1]
        lat0 <- PPole[1,2]
      }
      if(plot==T){
        if(on_plot==FALSE){
          plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
               xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
          sph_ortho(lat=lat0,long=lon0,grid=grid,coast=coast, title=title)
        }
      }

      coord <- as.data.frame(lon0)
      coord$lat0 <- lat0
      VGP$x <- c2x(VGP$lon,VGP$lat)
      VGP$y <- c2y(VGP$lon,VGP$lat)
      VGP$cut <- cut(VGP$lon,VGP$lat)

      #select symbol
      if(symbol=="c") {pch <- 21}
      else if(symbol=="s") {pch <- 22}
      else if(symbol=="d") {pch <- 23}
      else if(symbol=="t") {pch <- 24}

      if(plot==T){
        #plot points if requested
        if(plot_vgp==T){
          points(VGP$x,VGP$y,pch=pch,
                 col=col_sym_out,
                 bg=ifelse(VGP$cut>0,col,"white"), cex=cex)
        }
      }

      #uses PmagDiR for plotting A95
      if(plot==T){
        if(A95==TRUE){
          plot_PA95(lon = PPole[1,1],lat = PPole[1,2],A = PPole[1,3],lon0 = lon0,lat0 = lat0,on_plot = TRUE,symbol = symbol,col_f = col_f,size=1.5)
        }
      }

      #calculate bootstreapped VGPs
      if(B95==TRUE){
        bootlonlat <- as.data.frame(matrix(ncol = 2,nrow = 0))
        n <- 0
        #number of bootrstrap
        nb <- nb
        repeat{
          n <- n+1
          VGPb <- boots_DI(VGP)
          VGPb_av <- fisher(VGPb)
          blon <- VGPb_av[1,1]
          blat <- VGPb_av[1,2]
          blonlat <- as.data.frame(t(c(blon,blat)))
          bootlonlat <- rbind(bootlonlat,blonlat)
          x <- c2x(blon,blat)
          y <- c2y(blon,blat)
          cutt <- cut(blon,blat)

          #update progress bar of VGP calculated from directions internally
          if(VGPint==1){
            updateProgressBar(
              id="vgpboot",
              value=n,total=nb,
            )
          }else if(VGPint==2){
            #Update progress bar of VGP external
            updateProgressBar(
              id="Ext_vgpboot",
              value=n,total=nb,
            )
          }else if(VGPint==3){
            updateProgressBar(
              id="Mvgpboot",
              value=n,total=nb,
            )
          }else if(VGPint==4){
            updateProgressBar(
              id="SVGPboot",
              value=n,total=nb,
            )
          }

          if(n>=nb) break
          #ADD bootstrapping counter
        }
        colnames(bootlonlat) <- c("vgp_lon","vgp_lat")
        #calculate angluar distances
        bootlonlat$Plon <- rep(PPole[1,1])
        bootlonlat$Plat <- rep(PPole[1,2])
        bootlonlat$delta <- abs(bootlonlat$vgp_lon-bootlonlat$Plon)
        bootlonlat$diff <- r2d(acos((sin(d2r(bootlonlat$vgp_lat))*sin(d2r(bootlonlat$Plat)))+
                                      (cos(d2r(bootlonlat$vgp_lat))*cos(d2r(bootlonlat$Plat))*cos(d2r(bootlonlat$delta)))))
        ang_dis <- as.data.frame(bootlonlat$diff)
        ang_dis <- (ang_dis[order(ang_dis[,1]),])
        conf <- 0.95
        Uconf <- round(nb*conf,digits=0)
        angular_conf <- ang_dis[Uconf]

        #add cartesian coordinates
        bootlonlat$x <- c2x(bootlonlat$vgp_lon,bootlonlat$vgp_lat)
        bootlonlat$y <- c2y(bootlonlat$vgp_lon,bootlonlat$vgp_lat)
        bootlonlat$cutt <- cut(bootlonlat$vgp_lon,bootlonlat$vgp_lat)

        if(plot==T){
          #plot bootstrapped data
          points(bootlonlat$x,bootlonlat$y,pch=ifelse(bootlonlat$cutt>0,16,1),col=col_boot)
          PPole_x <- c2x(PPole[1,1],PPole[1,2])
          PPole_y <- c2y(PPole[1,1],PPole[1,2])
          PPole_cut <- cut(PPole[1,1],PPole[1,2])
          #plot average
          if(PPole_cut>0){
            points(PPole_x,PPole_y, pch=pch,cex=1.5, col="black",
                   bg= col_f)
          }else{
            points(PPole_x,PPole_y, pch=pch,cex=1.5, col="black",
                   bg= "white")
          }
        }
      }

      Plot_VGP_result <- list(0)
      Paleopole_F <- cbind(PPole[,4],PPole[,1:3],PPole[,6])
      colnames(Paleopole_F) <- c("N","Long","Lat","A95","K")
      rownames(Paleopole_F) <- input$fileN_VGP
      if(B95==TRUE){
        Paleopole_B <- data.frame(matrix(ncol=4,nrow=1))
        Paleopole_B[1,1] <- PPole[1,4]
        Paleopole_B[1,2:3] <- PPole[1,1:2]
        Paleopole_B[1,4] <- angular_conf
        colnames(Paleopole_B) <- c("N","Long","Lat","B95")
        rownames(Paleopole_B) <- input$fileN_VGP
      }

      Plot_VGP_result[[1]] <- Paleopole_F
      if(B95==TRUE){
        Plot_VGP_result[[2]] <- Paleopole_B
        Plot_VGP_result[[3]] <- bootlonlat[,1:2]
      }
      return(Plot_VGP_result)
    }

    #create reactive value
    Pole <- reactiveValues(FishPole = NULL)
    VGP <- reactiveValues(VGP_list=NULL)

    ### VGP calculation Part
    output$VGPplot <- renderPlot({
      #VGP_plot()
      DI<- fix_DI(input_file())
      if(input$dirs_vgp==1){DIrs <- DI}
      else if(input$dirs_vgp==2){               #save EI directions as reactive file and use it here!!!
        DIrs <- ffind_boot_S(DI,bootstrap = 1)
        DIrs <- DIrs[[1]]
      }
      #create file with all VGPs
      reactiveValuesToList(VGP)
      #VGP_list <- list()
      VGP$VGP_list[[1]] <- PmagDiR::VGP_DI(DI = DIrs,export = FALSE,Prnt = FALSE,lat = input$lat, long = input$long,type="VGPsN")
      VGP$VGP_list[[2]] <- PmagDiR::VGP_DI(DI = DIrs,export = FALSE,Prnt = FALSE,lat = input$lat, long = input$long,type="VGPs")
      VGP$VGP_list[[3]] <- PmagDiR::VGP_DI(DI = DIrs,export = FALSE,Prnt = FALSE,lat = input$lat, long = input$long,type="VGPsR")
      VGPS <- VGP$VGP_list[[1]]

      #select color
      if(input$vgpscolor==1) vgpscolor <- "black"
      if(input$vgpscolor==2) vgpscolor <- "blue"
      if(input$vgpscolor==3) vgpscolor <- "green"
      if(input$vgpscolor==4) vgpscolor <- "pink"
      if(input$vgpscolor==5) vgpscolor <- "purple"
      if(input$vgpscolor==6) vgpscolor <- "brown"
      if(input$vgpscolor==7) vgpscolor <- "red"
      if(input$vgpscolor==8) vgpscolor <- "yellow"
      if(input$vgpscolor==9) vgpscolor <- "cyan"
      if(input$vgpscolor==10) vgpscolor <- "gray"
      if(input$vgpscolor==11) vgpscolor <- "white"
      assign("VGP_color", vgpscolor, envir = MVGP_temp)

      #select symbol
      if(input$vgpssymbol==1) vgpssymbol <- "c"
      if(input$vgpssymbol==2) vgpssymbol <- "s"
      if(input$vgpssymbol==3) vgpssymbol <- "d"
      if(input$vgpssymbol==4) vgpssymbol <- "t"
      assign("VGP_symbol", vgpssymbol, envir = MVGP_temp)

      #define name with locality
      vgp_temp_stat_name <- paste(dirsFileName,"stat_temp",nrow(VGPS),input$lat,input$long,input$coord,input$vgpbootn,input$dirs_vgp,input$cutoff,sep = "_")
      if(input$plotA95==3 && exists(vgp_temp_stat_name,envir = .GlobalEnv)==FALSE) {
        boot_run <- TRUE
      }else{boot_run <- FALSE}
      #plot data and save stat
      Pole$FishPole <- plot_VGP_S(VGP = VGPS,
                                  lat = input$centerlat,
                                  long = input$centerlong,
                                  auto_cent = ifelse(input$centercoord==1,TRUE,FALSE),
                                  coast = ifelse(input$coastyesno==1,TRUE,FALSE),
                                  col = vgpscolor,symbol = vgpssymbol,nb=input$vgpbootn,
                                  A95 = ifelse(input$plotA95==2,TRUE,FALSE),
                                  B95 = boot_run,
                                  VGPint = 1)
      if(boot_run==TRUE){assign(x = vgp_temp_stat_name,value = Pole$FishPole,envir = .GlobalEnv)}

      #if bootstrap is selected but already exist, plot bootstrapped data
      if(input$plotA95==3 && exists(vgp_temp_stat_name,envir = .GlobalEnv)==TRUE) {
        plot_VGP_S(VGP = .GlobalEnv[[vgp_temp_stat_name]][[3]][,1:2],
                   lat = input$centerlat,
                   long = input$centerlong,
                   auto_cent = ifelse(input$centercoord==1,TRUE,FALSE),
                   on_plot = TRUE, col = rgb(1,0,0,0.15),col_sym_out = rgb(1,0,0,0.15))
        PmagDiR::plot_PA95(lon = .GlobalEnv[[vgp_temp_stat_name]][[2]][1,2],
                           lat = .GlobalEnv[[vgp_temp_stat_name]][[2]][1,3],
                           A = 0,
                           lat0 = ifelse(input$centercoord==2,input$centerlat,.GlobalEnv[[vgp_temp_stat_name]][[2]][1,3]),
                           lon0 = ifelse(input$centercoord==2,input$centerlong,.GlobalEnv[[vgp_temp_stat_name]][[2]][1,2]),
                           col_f = "cyan",
                           symbol = vgpssymbol,
                           size = 1.5,on_plot = T)
      }

      #make sure the pole reactive file is complete
      if(input$plotA95==3){Pole$FishPole <- .GlobalEnv[[vgp_temp_stat_name]]}

      output$fishpole <- renderTable({
        if(input$plotA95==1) {Polestat <- NULL}
        if(input$plotA95==2) {Polestat <- Pole$FishPole[[1]]}
        if(input$plotA95==3) {Polestat <- Pole$FishPole[[2]]}
        Polestat
      }, digits = 1,align = "c", rownames = T)
      #warning for geographic coordinates
      geowarn <- ifelse(input$coord==1,"WARNING: data are in Geographic coordinates","")
      output$geowarning <- renderText({geowarn})

      #warning for site coordinates not set
      coordwarn <- ifelse(input$lat==0 | input$long==0, "WARNING: site latitude and longitude are set as zero","")
      output$coordwarning <- renderText({coordwarn})

      VGPs_plot <- recordPlot()
      #Export graphic
      output$VGPs_G <- downloadHandler(
        filename = function() {
          paste(input$fileN_VGP,"_VGPs_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 11,height = 11)
          replayPlot(VGPs_plot)
          dev.off()
        }
      )
      output$VGPs_S <- downloadHandler(
        filename = function() {
          paste(input$fileN_VGP,"_POLE_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          if(input$plotA95==1){write.csv(NULL,file)}
          else if(input$plotA95==2){
            rownames(Pole$FishPole[[1]]) <- input$fileN_VGP
            write.csv(round(Pole$FishPole[[1]], digits = 2),file)
          }
          else if(input$plotA95==3){
            bootPole <- .GlobalEnv[[vgp_temp_stat_name]][[2]]
            rownames(bootPole) <- input$fileN_VGP
            write.csv(round(bootPole, digits = 2),file)
          }
        }
      )
      output$VGPs_Exp <- downloadHandler(
        filename = function(){
          paste(input$fileN_VGP,
                if(input$VGPtype==1){"_VGPs_SingleMod_"}
                else if(input$VGPtype==2){"_VGP_Rev_"}
                else if(input$VGPtype==3){"_VGP_Rot_"},
                Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          if(input$VGPtype==1){write.csv(round(VGP$VGP_list[[1]], digits = 2),file, row.names = F)}
          else if(input$VGPtype==2){write.csv(round(VGP$VGP_list[[2]], digits = 2),file, row.names = F)}
          else if(input$VGPtype==3){write.csv(round(VGP$VGP_list[[3]], digits = 2),file, row.names = F)}
        }
      )
    }, width = 700, height = 700)

    #send current VGP to list and reactive VGP_saved_list
    observeEvent(input$saveVGP,{
      req(VGP$VGP_list[[1]])
      VGP <- VGP$VGP_list[[1]]
      colnames(VGP) <- c("Long","Lat")
      VGP$color <- MVGP_temp[["VGP_color"]]
      VGP$symbol <- MVGP_temp[["VGP_symbol"]]

      VGPsaved <- list()
      VGPsaved[[1]] <- VGP
      #if statistics does not exists, it makes it
      if(length(Pole$FishPole)<=1){
        VGPsaved[[2]] <- plot_VGP_S(VGP = VGP,A95 = T,B95 = T,nb = input$vgpbootn,VGPint = 1,plot=F)
      }else{VGPsaved[[2]] <- Pole$FishPole}
      VGP_saved$list[[input$fileN_VGP]] <- VGPsaved

      MVGP_list$vgps <- sites_list()
    })

    #preliminary table for loaded VGP
    output$Ext_VGP_list3 <- renderTable({
      if(length(VGP_saved$list)==0){
        MVGP_list$vgps <- NULL
      }else{
        MVGP_list$vgps <- sites_list()
      }
      MVGP_list$vgps
    }, digits = 1,align = "l", rownames = T)

    ### External and Multiple vgps part
    #creates reactive value for checking if file is uploaded
    values <- reactiveValues(mvgp = NULL)
    #check for uploaded file
    observeEvent(input$vgpfile,{values$mvgp <- "uploaded"})
    #reset upload if requested
    observeEvent(input$resetMVGP,{values$mvgp <- "reset"})
    #read file if present, reset if requested
    VGPfile <- reactive({
      if (is.null(values$mvgp)) {
        return(NULL)
      } else if (values$mvgp == 'uploaded') {
        read.csv(file = input$vgpfile$datapath)
      } else if (values$mvgp == 'reset') {
        return(NULL)
      }
    })

    #function that apply cutoff to VGP
    cut_VGP <- function(VGP){
      #degree to radians and vice versa
      d2r <- function(x) {x*(pi/180)}
      r2d <- function(x) {x*(180/pi)}
      colnames(VGP) <- c("Long","Lat")
      N <- nrow(VGP)
      #calculate average
      Pole <- PmagDiR::fisher(VGP)
      VGP$Plong <- rep(Pole[1,1])
      VGP$Plat <- rep(Pole[1,2])
      #calculate distance pole - VGP
      VGP$delta <- abs(VGP$Long-VGP$Plong)
      VGP$diff <- r2d(acos((sin(d2r(VGP$Lat))*sin(d2r(VGP$Plat)))+
                             (cos(d2r(VGP$Lat))*cos(d2r(VGP$Plat))*cos(d2r(VGP$delta)))))

      #apply cutoff
      if(input$VGP_ext_cut_type==2){
        #vandamme filtering calculation
        ASD <- sqrt(sum(((VGP$diff)^2)/(N-1)))
        A <- (1.8*ASD)+5
      }else if(input$VGP_ext_cut_type==3){
        A <- input$VGP_ext_cutoff
      }
      VGPcut <- as.numeric(which(VGP$diff>A), arr.ind = TRUE)
      if(length(VGPcut)!=0) {VGP <- VGP[-VGPcut,]}
      VGP <- VGP[,1:2]
      return(VGP)
    }

    #create service environment
    MVGP_temp <- new.env()
    assign("MVGP_temp", MVGP_temp, envir = .GlobalEnv)

    #plot current vgp and save plotting detail in temporary environment
    plot_current_VGP <- function(){
      req(VGPfile())
      VGP <- VGPfile()
      colnames(VGP) <- c("Long","Lat")
      #apply cutoff if requested
      if(input$VGP_ext_cut_type!=1){VGP <- cut_VGP(VGP = VGP)}

      #choose color
      if(input$MVGPcolor==1) MVGPcolor <- "black"
      if(input$MVGPcolor==2) MVGPcolor <- "blue"
      if(input$MVGPcolor==3) MVGPcolor <- "green"
      if(input$MVGPcolor==4) MVGPcolor <- "pink"
      if(input$MVGPcolor==5) MVGPcolor <- "purple"
      if(input$MVGPcolor==6) MVGPcolor <- "brown"
      if(input$MVGPcolor==7) MVGPcolor <- "red"
      if(input$MVGPcolor==8) MVGPcolor <- "yellow"
      if(input$MVGPcolor==9) MVGPcolor <- "cyan"
      if(input$MVGPcolor==10) MVGPcolor <- "gray"
      if(input$MVGPcolor==11) MVGPcolor <- "white"


      #send color to temporary environment
      assign("MVGPcolor", MVGPcolor, envir = MVGP_temp)      #eliminate environment and use only reactive values

      #choose symbol
      if(input$MVGPsymbol==1) MVGPsymbol <- "c"
      if(input$MVGPsymbol==2) MVGPsymbol <- "s"
      if(input$MVGPsymbol==3) MVGPsymbol <- "d"
      if(input$MVGPsymbol==4) MVGPsymbol <- "t"
      #send symbol to temporary environment
      assign("MVGPsymbol", MVGPsymbol, envir = MVGP_temp)


      #plot vgp and save statistic on temporary file
      #ask if bootrstrapped stat with same localitiy name an VGPs number already exists and does not perform it iyes
      temp_stat_name <- paste(input$VGP_ext_sitename,"_stat_temp_",nrow(VGP),"_nb_",input$MVGPnb_ext,sep = "")
      if(input$MVGP_stat==3 && exists(temp_stat_name,envir = MVGP_temp)==FALSE) {
        boot_run <- TRUE
      }else{boot_run <- FALSE}
      MVGP_Pole <- plot_VGP_S(VGP = VGP,lat = input$VGP_ext_clat,
                              long = input$VGP_ext_clong,
                              auto_cent = ifelse(input$VGP_ext_center==1,TRUE,FALSE),
                              coast = ifelse(input$VGP_ext_coast==1,TRUE,FALSE),
                              col = MVGPcolor, symbol = MVGPsymbol,nb = input$MVGPnb,
                              A95 = ifelse(input$MVGP_stat==2,TRUE,FALSE),
                              B95 = boot_run,
                              VGPint=2)

      #if bootstrap is selected but already exist, plot bootstrapped data
      if(input$MVGP_stat==3 && exists(temp_stat_name,envir = MVGP_temp)==TRUE) {
        plot_VGP_S(VGP = MVGP_temp[[temp_stat_name]][[3]][,1:2],
                   lat = input$VGP_ext_clat,
                   long = input$VGP_ext_clong,
                   auto_cent = ifelse(input$VGP_ext_center==1,TRUE,FALSE),
                   on_plot = TRUE, col = rgb(1,0,0,0.15),col_sym_out = rgb(1,0,0,0.15))
        plot_VGP_S(VGP = MVGP_temp[[temp_stat_name]][[2]][1,2:3],
                   lat = input$VGP_ext_clat,
                   long = input$VGP_ext_clong,
                   auto_cent = ifelse(input$VGP_ext_center==1,TRUE,FALSE),
                   on_plot = TRUE, cex=1.5, col= "cyan",symbol = "d")
      }
      #save statistic in a file with name of locality and number of VGPs, so bootstrap is not performed again if already exists
      if(input$MVGP_stat==3){
        if(exists(temp_stat_name,envir = MVGP_temp)==FALSE){assign(x = temp_stat_name,value = MVGP_Pole,envir = MVGP_temp)}
      }
      final_pole <- MVGP_temp[[temp_stat_name]]
      #send temporary file to global environment
      assign("MVGP_Pole",final_pole, envir = MVGP_temp)

      #populate and plot CURRENT VGP pole statistic
      output$MVGPpolestat <- renderTable({
        if(input$MVGP_stat==1) {MVGPpolestat <- NULL}
        if(input$MVGP_stat==2) {
          MVGPpolestat <- MVGP_Pole[[1]]
          rownames(MVGPpolestat) <- input$VGP_ext_sitename
        }
        if(input$MVGP_stat==3) {
          MVGPpolestat <- MVGP_temp[[temp_stat_name]][[2]]
          rownames(MVGPpolestat) <- input$VGP_ext_sitename
        }
        if(values$mvgp == 'reset'){MVGPpolestat <- NULL}
        MVGPpolestat
      }, digits = 1,align = "l", rownames = T)

      #export locality level stat
      output$VGP_site_stat <- downloadHandler(
        filename = function() {
          paste(input$VGP_ext_sitename,"_POLE_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          if(input$MVGP_stat==1){write.csv(NULL,file)}
          else if(input$MVGP_stat==2){
            POLE <- MVGP_Pole[[1]]
            rownames(POLE) <- paste(input$VGP_ext_sitename,"_fisher",sep="")
            write.csv(round(POLE, digits = 2),file)
          }
          else if(input$MVGP_stat==3){
            POLE <- MVGP_temp[[temp_stat_name]][[2]]
            rownames(POLE) <- paste(input$VGP_ext_sitename,"_bootstrap",sep="")
            write.csv(round(POLE, digits = 2),file)
          }
        }
      )
    }

    #function that built list of sites
    sites_list <- function(){
      sites <- ls(VGP_saved$list)
      MVGP_list_t <- data.frame(matrix(nrow=0,ncol = 8))
      colnames(MVGP_list_t) <- c("Loc","Col", "Sym","N","Lon","Lat","A95","B95")
      #populate list
      for(i in 1:length(sites)){
        newTabLine <- data.frame(t(c(sites[i],
                                     VGP_saved$list[[sites[i]]][[1]][1,3],
                                     VGP_saved$list[[sites[i]]][[1]][1,4],
                                     nrow(VGP_saved$list[[sites[i]]][[1]]),
                                     round(VGP_saved$list[[sites[i]]][[2]][[1]][1,2], digits=1),
                                     round(VGP_saved$list[[sites[i]]][[2]][[1]][1,3], digits=1),
                                     round(VGP_saved$list[[sites[i]]][[2]][[1]][1,4], digits=1),
                                     round(VGP_saved$list[[sites[i]]][[2]][[2]][1,4], digits=1))))
        colnames(newTabLine) <- c("Loc.","Col", "Sym","N","Lon","Lat","A95","B95")
        MVGP_list_t <- rbind(MVGP_list_t,newTabLine)
      }
      return(MVGP_list_t)
    }

    #create reactive files for saving external VGPs
    VGP_saved <- reactiveValues(list = NULL)
    MVGP_list <- reactiveValues(vgps = NULL)

    #send current VGP to list and reactive VGP_saved_list
    observeEvent(input$saveMVGP,{
      req(VGPfile())
      VGP <- VGPfile()
      colnames(VGP) <- c("Long","Lat")
      #apply cutoff if requested
      if(input$VGP_ext_cut_type!=1){VGP <- cut_VGP(VGP = VGP)}
      VGP$color <- MVGP_temp[["MVGPcolor"]]
      VGP$symbol <- MVGP_temp[["MVGPsymbol"]]

      VGPsaved <- list()
      VGPsaved[[1]] <- VGP
      #if statistics does not exists, it makes it
      if(length(MVGP_temp[["MVGP_Pole"]])<=1){
        VGPsaved[[2]] <- plot_VGP_S(VGP = VGP,A95 = T,B95 = T,nb = input$MVGPnb_ext,VGPint = 2,plot=F)
      }else{VGPsaved[[2]] <- MVGP_temp[["MVGP_Pole"]]}
      VGP_saved$list[[input$VGP_ext_sitename]] <- VGPsaved

      MVGP_list$vgps <- sites_list()
    })

    #delete VGP from reactive lists
    observeEvent(input$deletevgp,{
      e <- input$MVGPlist_rows_selected
      for(i in e){
        #eliminate selected names from list
        VGP_saved$list <- VGP_saved$list[names(VGP_saved$list) %in% MVGP_list$vgps[i,1]==F]
      }
      if(length(VGP_saved$list)==0){
        MVGP_list$vgps <- NULL
      }else{
        MVGP_list$vgps <- sites_list()
      }
    })

    #send figures to ui external VGP
    output$Ext_VGP_plot <- renderPlot({
      plot_current_VGP()

      #record plot
      Ext_VGP_plot <- recordPlot()

      #Export graphic
      output$Ext_VGP_G <- downloadHandler(
        filename = function() {
          paste(input$VGP_ext_sitename,"_MultiVGP_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 11,height = 11)
          replayPlot(Ext_VGP_plot)
          dev.off()
        }
      )},width = 700, height = 700)

    #preliminary table for loaded VGP
    output$Ext_VGP_list <- renderTable({
      if(length(VGP_saved$list)==0){
        MVGP_list$vgps <- NULL
      }else{
        MVGP_list$vgps <- sites_list()
      }
      MVGP_list$vgps
    }, digits = 1,align = "l", rownames = T)

    #interactive table for loaded VGP
    output$MVGPlist <- DT::renderDataTable(MVGP_list$vgps, server = F)

    #interactive table for loaded VGP in Euler rotation page
    output$MVGPlist2 <- DT::renderDataTable(MVGP_list$vgps, server = F,selection="single")

    #plot VGPs selected from list for euler rotation
    plot_selected_VGP_eul <- function(r,lat0,lon0){
      #selects what to plot
      if(input$eulPlotType==1){
        plot_VGP_S(VGP=VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][,1:2],
                   lat = lat0,
                   long = lon0,
                   auto_cent = F, on_plot = T,
                   col = VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][1,3],
                   symbol = VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][1,4])
      }else if(input$eulPlotType==2){
        #plot fisher means and not VGPs if requested
        PmagDiR::plot_PA95(lon = VGP_saved$list[[MVGP_list$vgps[r,1]]][[2]][[1]][1,2],
                           lat = VGP_saved$list[[MVGP_list$vgps[r,1]]][[2]][[1]][1,3],
                           A = VGP_saved$list[[MVGP_list$vgps[r,1]]][[2]][[1]][1,4],
                           lat0 = lat0,
                           lon0 = lon0,
                           col_A = rgb(0,0,1,0.15),
                           symbol = VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][1,4],
                           col_f = VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][1,3],
                           on_plot = T)
      }else if(input$eulPlotType==3){
        plot_VGP_S(VGP=VGP_saved$list[[MVGP_list$vgps[r,1]]][[2]][[3]][,1:2],
                   lat = lat0,
                   long = lon0,
                   auto_cent = F, on_plot = T,
                   col = rgb(1,0,0,0.1),col_sym_out = rgb(1,0,0,0.1),
                   symbol = "c")
        PmagDiR::plot_PA95(lon = VGP_saved$list[[MVGP_list$vgps[r,1]]][[2]][[1]][1,2],
                           lat = VGP_saved$list[[MVGP_list$vgps[r,1]]][[2]][[1]][1,3],
                           A = 0,
                           lon0 = lon0,lat0 = lat0,
                           col_f = VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][1,3],
                           symbol = VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][1,4],
                           size = 1.2,on_plot = T)
      }
    }

    #creates reactive file
    vgprotate <- reactiveValues(new=NULL)

    #function that rotate the data
    observeEvent(input$eulerrot,{
      r <- input$MVGPlist2_rows_selected
      if(length(r)){
        vgprotate$new <- PmagDiR::rot_DI(Lonlat = VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]][,1:2],
                                         P_long = input$eul_long,P_lat = input$eul_lat,rot = input$eul_rot)
      }
    })

    #creates reactive file with no values to avoid problems in plotting module
    vgprotate <- reactiveValues(new=NULL)

    #function rotating the data and stats
    observeEvent(input$eulerrot,{
      r <- input$MVGPlist2_rows_selected
      if(length(r)){
        #choose color of average
        if(input$eulcolor==1) eulcolor <- "black"
        if(input$eulcolor==2) eulcolor <- "blue"
        if(input$eulcolor==3) eulcolor <- "green"
        if(input$eulcolor==4) eulcolor <- "pink"
        if(input$eulcolor==5) eulcolor <- "purple"
        if(input$eulcolor==6) eulcolor <- "brown"
        if(input$eulcolor==7) eulcolor <- "red"
        if(input$eulcolor==8) eulcolor <- "yellow"
        if(input$eulcolor==9) eulcolor <- "cyan"
        if(input$eulcolor==10) eulcolor <- "gray"
        if(input$eulcolor==11) eulcolor <- "white"

        #select symbol of average
        if(input$eulsymbol==1) eulsymbol <- "c"
        if(input$eulsymbol==2) eulsymbol <- "s"
        if(input$eulsymbol==3) eulsymbol <- "d"
        if(input$eulsymbol==4) eulsymbol <- "t"

        vgp_rot <- list()

        #copy file to rotate in new list
        vgp_rot[[1]] <- VGP_saved$list[[MVGP_list$vgps[r,1]]][[1]]
        #rotate VGPs
        vgp_rot[[1]][,1:2] <- PmagDiR::rot_DI(Lonlat = vgp_rot[[1]][,1:2],
                                              P_long = input$eul_long,P_lat = input$eul_lat,rot = input$eul_rot)
        #copy color and symbol to file for saving it in list
        vgp_rot[[1]][,3] <- rep(eulcolor)
        vgp_rot[[1]][,4] <- rep(eulsymbol)


        #copy original VGPs stats to be rotated
        vgp_rot[[2]] <- VGP_saved$list[[MVGP_list$vgps[r,1]]][[2]]
        #rotate fisher
        vgp_rot[[2]][[1]][1,2:3] <- PmagDiR::rot_DI(Lonlat = vgp_rot[[2]][[1]][1,2:3],
                                                    P_long = input$eul_long,P_lat = input$eul_lat,rot = input$eul_rot)
        #rotate B95
        vgp_rot[[2]][[2]][1,2:3] <- PmagDiR::rot_DI(Lonlat = vgp_rot[[2]][[2]][1,2:3],
                                                    P_long = input$eul_long,P_lat = input$eul_lat,rot = input$eul_rot)
        #rotate bootstreapped dirs
        vgp_rot[[2]][[3]][,1:2] <- PmagDiR::rot_DI(Lonlat = vgp_rot[[2]][[3]][,1:2],
                                                   P_long = input$eul_long,P_lat = input$eul_lat,rot = input$eul_rot)
        #transfer data to reactive file
        vgprotate$new <- vgp_rot
      }
    })

    #function adding rotated data to list
    observeEvent(input$eulersave,{
      VGP_saved$list[[input$eul_name]] <- vgprotate$new
      sites <- ls(VGP_saved$list)
      MVGP_list$vgps <- sites_list()
    })

    #function deleting rotated data
    observeEvent(input$eulerdel,{
      vgprotate$new <- NULL
    })

    #build rotated stat temporary table
    output$eul_temp_table <- renderTable({
      #check if rotated file exists
      if(is.null(vgprotate$new)==F){
        rot_temp_stat <- data.frame(matrix(nrow=0,ncol = 8))
        colnames(rot_temp_stat) <- c("Loc","Col", "Sym","N","Lon","Lat","A95","B95")
        rot_temp_stat[1,1] <- input$eul_name
        rot_temp_stat[1,2] <- vgprotate$new[[1]][1,3]
        rot_temp_stat[1,3] <- vgprotate$new[[1]][1,4]
        rot_temp_stat[1,4] <- nrow(vgprotate$new[[1]])
        rot_temp_stat[1,5:6] <- vgprotate$new[[2]][[1]][1,2:3]
        rot_temp_stat[1,7] <- vgprotate$new[[2]][[1]][1,4]
        rot_temp_stat[1,8] <- vgprotate$new[[2]][[2]][1,4]
        rot_temp_stat
      }
    },digits = 1)

    #send figure to ui
    output$eulerplot <- renderPlot({
      r <- input$MVGPlist2_rows_selected

      if(length(r)>1){r <- r[1]}
      if(length(r)){
        centerLat = ifelse(input$eul_center==1, VGP_saved$list[[MVGP_list$vgps[r[1],1]]][[2]][[1]][1,3],input$eul_clat)
        centerLong = ifelse(input$eul_center==1, VGP_saved$list[[MVGP_list$vgps[r[1],1]]][[2]][[1]][1,2],input$eul_clong)
      }else{
        centerLat <- input$eul_clat
        centerLong <- input$eul_clong
      }
      if(length(r)){

        #plot empty spherical orthographic
        PmagDiR::sph_ortho(lat = centerLat,long = centerLong,
                           coast = ifelse(input$eul_coast==1,TRUE,FALSE))

        plot_selected_VGP_eul(r = r,lat0 = centerLat,lon0 = centerLong)
        if(is.null(vgprotate$new)==F){
          if(input$eulPlotType==1){
            plot_VGP_S(VGP = vgprotate$new[[1]],
                       lat = centerLat,
                       long = centerLong,
                       col= vgprotate$new[[1]][1,3],
                       symbol = vgprotate$new[[1]][1,4],
                       on_plot = T, auto_cent = F)
          }else if(input$eulPlotType==2){
            PmagDiR::plot_PA95(lon = vgprotate$new[[2]][[1]][1,2],
                               lat = vgprotate$new[[2]][[1]][1,3],
                               A = vgprotate$new[[2]][[1]][1,4],
                               lat0 = centerLat,
                               lon0 = centerLong,
                               col_A = rgb(1,0,1,0.15),
                               symbol = vgprotate$new[[1]][1,4],
                               col_f = vgprotate$new[[1]][1,3],
                               on_plot = T)
          }else if(input$eulPlotType==3){
            plot_VGP_S(VGP=vgprotate$new[[2]][[3]][,1:2],
                       lat = centerLat,
                       long = centerLong,
                       auto_cent = F, on_plot = T,
                       col = rgb(1,0,0,0.1),col_sym_out = rgb(1,0,0,0.1),
                       symbol = "c")
            PmagDiR::plot_PA95(lon = vgprotate$new[[2]][[2]][1,2],
                               lat = vgprotate$new[[2]][[2]][1,3],
                               A = 0,
                               lon0 = centerLong,
                               lat0 = centerLat,
                               col_f = vgprotate$new[[1]][1,3],
                               symbol = vgprotate$new[[1]][1,4],
                               size = 1.2,on_plot = T)
          }
        }
      }
      #record plot
      euler_plot <- recordPlot()

      #Export graphic
      output$euler_G <- downloadHandler(
        filename = function() {
          paste(input$eul_name,"_rot_vgp_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 11,height = 11)
          replayPlot(euler_plot)
          dev.off()
        }
      )

      #export rotated vgps
      output$euler_vgp <- downloadHandler(
        filename = function() {
          paste(input$eul_name,"_rot_vgp_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          if(is.null(vgprotate$new)==F){
            write.csv(round(vgprotate$new[[1]][,1:2], digits = 2),file,row.names = F)
          }
        }
      )

    },width = 700, height = 700)

    ######## SEND PLOT TO UI OF BOTH MULTIPLE VGP PARTS #########
    #plot VGPs selected from list
    plot_selected_VGP <- function(s,lat0,lon0){

      #set bootfile name
      for(i in s){
        #selects what to plot
        if(input$MVGPsPlotType==1){
          plot_VGP_S(VGP=VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][,1:2],
                     lat = lat0,
                     long = lon0,
                     auto_cent = F, on_plot = T,
                     col = VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][1,3],
                     symbol = VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][1,4])
        }else if(input$MVGPsPlotType==2){
          #plot fisher means and not VGPs if requested
          PmagDiR::plot_PA95(lon = VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][1,2],
                             lat = VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][1,3],
                             A = VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][1,4],
                             lat0 = lat0,
                             lon0 = lon0,
                             col_A = rgb(0,0,1,0.15),
                             symbol = VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][1,4],
                             col_f = VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][1,3],
                             on_plot = T)
        }else if(input$MVGPsPlotType==3){
          plot_VGP_S(VGP=VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[3]][,1:2],
                     lat = lat0,
                     long = lon0,
                     auto_cent = F, on_plot = T,
                     col = rgb(1,0,0,0.1),col_sym_out = rgb(1,0,0,0.1),
                     symbol = "c")
          PmagDiR::plot_PA95(lon = VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][1,2],
                             lat = VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][1,3],
                             A = 0,
                             lon0 = lon0,lat0 = lat0,
                             col_f = VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][1,3],
                             symbol = VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][1,4],
                             size = 1.2,on_plot = T)
        }

        #plot pole names if required
        if(input$MVGP_names_YN==2){
          #functions converting degree and radians
          d2r <- function(x) {x*(pi/180)}
          r2d <- function(x) {x*(180/pi)}
          #functions converting long & lat to xy
          c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
          c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
          #cut is cosin of c, when negative is behind projections, needs to be cut
          cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}

          #define coordinate for name of pole
          x <- c2x(VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][[1,2]],VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][[1,3]])
          y <- c2y(VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][[1,2]],VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][[1,3]])
          name <- MVGP_list$vgps[i,1]
          #plot names
          text(x=x, y=y,pos=3,substitute(paste(bold(name))), cex= 1.2)
        }
      }
    }

    #function that creates plots for Multiple VGP analysis
    all_poles_plotter <- function(){
      s <- input$MVGPlist_rows_selected
      ext <- input$EP_list_rows_selected
      #if no pole is selected from list, read center from input always
      if(length(s)==0 && length(ext)==0){
        centerLat <- input$MVGP_clat
        centerLong <- input$MVGP_clong
      }else if(length(s)){
        centerLat = ifelse(input$MVGP_center==1, VGP_saved$list[[MVGP_list$vgps[s[1],1]]][[2]][[1]][1,3],input$MVGP_clat)
        centerLong = ifelse(input$MVGP_center==1, VGP_saved$list[[MVGP_list$vgps[s[1],1]]][[2]][[1]][1,2],input$MVGP_clong)
      } else if(length(ext)!=0 && length(s)==0){
        centerLat= ifelse(input$MVGP_center==1,Added_poles$list[ext[1],5] ,input$MVGP_clat)
        centerLong= ifelse(input$MVGP_center==1,Added_poles$list[ext[1],4] ,input$MVGP_clong)
      }
      #plot empty spherical orthographic
      PmagDiR::sph_ortho(lat = centerLat,long = centerLong,
                         coast = ifelse(input$MVGP_coast==1,TRUE,FALSE))

      #add APWP
      #read custom apwp file avoiding warning if not existing
      customAPWP <- reactive({
        #avoid warning if file is not loaded
        if (is.null(input$customAPWP)) {
          return(NULL)
        }
        # actually read the file
        read.csv(file = input$customAPWP$datapath)
      })

      #select APWP, either Vaes+2023 or Torsvik+2012
      if(input$APWP==2) {
        apwp <- "V23"
        frame <- as.numeric(input$frameV23)
      }
      if(input$APWP==3) {
        apwp <- "T12"
        frame <- as.numeric(input$frameT12)
      }

      #round to 10 if T12 or V23 are selected, to avoid problems in plotting pole using PmagDiR::plot_APWP
      if(input$APWP==2 || input$APWP==3){
        Y <- round(input$apwp_Y,-1)
        O <- round(input$apwp_O,-1)
      }else{
        Y <- input$apwp_Y
        O <- input$apwp_O
      }

      #plot apwp either Vaes+2023 or Torsvik+2012
      if(input$APWP==2 || input$APWP==3){
        PmagDiR::plot_APWP(APWP = apwp,
                           lat0 = centerLat,
                           lon0 = centerLong,
                           frame = frame,Shiny = T,Y = Y,O = O,size = 1.2,
                           coast = ifelse(input$MVGP_coast==1,TRUE,FALSE),Age_size = 1.5,on_plot = TRUE)
      }

      #plot custom APWP,
      if(input$APWP==4){
        req(customAPWP())
        c_apwp <- customAPWP()
        #select age interval depending on slide select
        c_apwp <- c_apwp[(c_apwp[,1]>=Y & c_apwp[,1]<=O),]
        #plot line connecting poles
        lin <- as.data.frame(c2x(c_apwp[,2],c_apwp[,3]))
        colnames(lin) <- "lx"
        lin$ly <- c2y(c_apwp[,2],c_apwp[,3])
        lines(lin$lx,lin$ly,cex=1)
        #plot poles
        for(i in 1:nrow(c_apwp)){
          PmagDiR::plot_PA95(lon = c_apwp[i,2],lat = c_apwp[i,3],A = c_apwp[i,4],
                             lon0 = centerLong,lat0 = centerLat,
                             col_f ="gray",col_A = rgb(1,0.9,0,0.30),on_plot = TRUE)
        }
        #plot min and max age
        text1 <- paste(c_apwp[1,1],"Ma")
        text2 <- paste(c_apwp[nrow(c_apwp),1], "Ma")
        text(x=lin[1,1], y=lin[1,2],pos=4,substitute(paste(bold(text1))), cex= 1.5)
        text(x=lin[nrow(lin),1], y=lin[nrow(lin),2],pos=4,substitute(paste(bold(text2))), cex= 1.5)
      }

      #functions converting degree and radians
      d2r <- function(x) {x*(pi/180)}
      r2d <- function(x) {x*(180/pi)}
      #functions converting long & lat to xy
      #convert to spherical coordinates using center Long and Lat as defined above
      c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-centerLong))}
      c2y <- function(lon,lat) {(cos(d2r(centerLat))*sin(d2r(lat)))-(sin(d2r(centerLat))*cos(d2r(lat))*cos(d2r(lon-centerLong)))}
      #cut is cosin of c, when negative is behind projections, needs to be cut
      cut <- function(lon,lat) {(sin(d2r(centerLat))*sin(d2r(lat)))+(cos(d2r(centerLat))*cos(d2r(lat))*cos(d2r(lon-centerLong)))}

      #if any vgp in table is selected it plots them, and the statistic, otherwise plots the currently loaded VGP
      if(length(s)){
        #choose color of average
        if(input$MVGP_aver_color==1) col_f <- "black"
        if(input$MVGP_aver_color==2) col_f <- "blue"
        if(input$MVGP_aver_color==3) col_f <- "green"
        if(input$MVGP_aver_color==4) col_f <- "pink"
        if(input$MVGP_aver_color==5) col_f <- "purple"
        if(input$MVGP_aver_color==6) col_f <- "brown"
        if(input$MVGP_aver_color==7) col_f <- "red"
        if(input$MVGP_aver_color==8) col_f <- "yellow"
        if(input$MVGP_aver_color==9) col_f <- "cyan"
        if(input$MVGP_aver_color==10) col_f <- "gray"
        if(input$MVGP_aver_color==11) col_f <- "white"

        #select symbol of average
        if(input$MVGP_aver_sym==1) MVGPaversym <- "c"
        if(input$MVGP_aver_sym==2) MVGPaversym <- "s"
        if(input$MVGP_aver_sym==3) MVGPaversym <- "d"
        if(input$MVGP_aver_sym==4) MVGPaversym <- "t"


        #no statistic
        if(input$MVGP_Pole_Stat==1){plot_selected_VGP(s=s,lat0 = centerLat,lon0 = centerLong)}

        #Fisher on Fisher algorithm
        if(input$MVGP_Pole_Stat==2){
          VGP_LonLat <- data.frame(matrix(ncol=2,nrow=0))
          colnames(VGP_LonLat) <- c("Long","Lat")
          for(i in s){
            LonLat_t <- VGP_saved$list[[MVGP_list$vgps[i,1]]][[2]][[1]][1,2:3]
            VGP_LonLat <- rbind(VGP_LonLat,LonLat_t)
          }
          MVGP_Fish_on_Fish <- fisher(VGP_LonLat)
          #rebuild results for table
          MVGP_ALL_Fisher <- data.frame(matrix(ncol=5, nrow = 1),row.names = input$fileN_MVGP)
          colnames(MVGP_ALL_Fisher) <- c("N","Long","Lat","A95","K")
          MVGP_ALL_Fisher[1,1] <- MVGP_Fish_on_Fish[1,4]
          MVGP_ALL_Fisher[1,2] <- MVGP_Fish_on_Fish[1,1]
          MVGP_ALL_Fisher[1,3] <- MVGP_Fish_on_Fish[1,2]
          MVGP_ALL_Fisher[1,4] <- MVGP_Fish_on_Fish[1,3]
          MVGP_ALL_Fisher[1,5] <- MVGP_Fish_on_Fish[1,6]

          plot_selected_VGP(s=s,lat0 = centerLat,lon0 = centerLong)
          PmagDiR::plot_PA95(lon = MVGP_Fish_on_Fish[1,1],lat = MVGP_Fish_on_Fish[1,2],A = MVGP_Fish_on_Fish[1,3],
                             lat0 = centerLat,
                             lon0 = centerLong,
                             on_plot = T,size = 1.5,col_f = col_f,symbol = MVGPaversym)
          #average pole name
          if(input$MVGP_names_YN==2){
            N_MVGP <- input$fileN_MVGP
            #add name
            text(x = c2x(MVGP_Fish_on_Fish[1,1],MVGP_Fish_on_Fish[1,2]),
                 y = c2y(MVGP_Fish_on_Fish[1,1],MVGP_Fish_on_Fish[1,2]),
                 pos=3,
                 substitute(paste(bold(N_MVGP))), cex= 1.5)
          }

          #fisher or bootstrap on VGP
        }else if(input$MVGP_Pole_Stat==3|input$MVGP_Pole_Stat==4){
          #calculate fisher or bootstrap of all VGPS
          VGP_LonLat <- data.frame(matrix(ncol=2,nrow=0))
          colnames(VGP_LonLat) <- c("Long","Lat")

          #set bootfile name for not repeating statistic
          MVGP_name <- "MVGP_stat"
          for(i in s){
            #fill the file name with all sites selected
            MVGP_name <- paste(MVGP_name,MVGP_list$vgps[i,1],sep="_")
            #combine all VGPs of selected sets
            LonLat_t <- VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]][1:nrow(VGP_saved$list[[MVGP_list$vgps[i,1]]][[1]]),1:2]
            VGP_LonLat <- rbind(VGP_LonLat,LonLat_t)
          }
          #calculate fisher
          if(input$MVGP_Pole_Stat==3){
            MVGP_Fish_on_VGP <- fisher(VGP_LonLat)
            #rebuild results for table
            MVGP_ALL_Fisher <- data.frame(matrix(ncol=5, nrow = 1),row.names = input$fileN_MVGP)
            colnames(MVGP_ALL_Fisher) <- c("N","Long","Lat","A95","K")
            MVGP_ALL_Fisher[1,1] <- MVGP_Fish_on_VGP[1,4]
            MVGP_ALL_Fisher[1,2] <- MVGP_Fish_on_VGP[1,1]
            MVGP_ALL_Fisher[1,3] <- MVGP_Fish_on_VGP[1,2]
            MVGP_ALL_Fisher[1,4] <- MVGP_Fish_on_VGP[1,3]
            MVGP_ALL_Fisher[1,5] <- MVGP_Fish_on_VGP[1,6]

            plot_selected_VGP(s=s,lat0 = centerLat,lon0 = centerLong)
            PmagDiR::plot_PA95(lon = MVGP_Fish_on_VGP[1,1],lat = MVGP_Fish_on_VGP[1,2],A = MVGP_Fish_on_VGP[1,3],
                               lat0 = centerLat,
                               lon0 = centerLong,
                               on_plot = T,size = 1.5,col_f = col_f, symbol = MVGPaversym)
            #add name if requested
            if(input$MVGP_names_YN==2){
              N_MVGP <- input$fileN_MVGP
              #add name
              text(x = c2x(MVGP_Fish_on_VGP[1,1],MVGP_Fish_on_VGP[1,2]),
                   y = c2y(MVGP_Fish_on_VGP[1,1],MVGP_Fish_on_VGP[1,2]),
                   pos=3,
                   substitute(paste(bold(N_MVGP))), cex= 1.5)
            }
          }
          #calculate bootstrap
          #build verifying name for existing bootstrap
          MVGP_temp_stat_name <- paste(MVGP_name,"_stat_temp_",nrow(VGP_LonLat),"_nb_",input$MVGPnb,sep = "")
          #if bootstrap is selected asks if file exists already in MVGP_temp envir
          if(input$MVGP_Pole_Stat==4 && exists(MVGP_temp_stat_name, envir = MVGP_temp)==FALSE){
            plot_selected_VGP(s=s,lat0 = centerLat,lon0 = centerLong)
            MVGP_boot_on_VGP <- plot_VGP_S(VGP = VGP_LonLat,
                                           lat = centerLat,
                                           lon = centerLong,
                                           plot_vgp = F,col_f = col_f,on_plot = T,B95 = T,nb =  input$MVGPnb,
                                           col_boot = rgb(0,0,1,0.15),auto_cent = F,symbol = MVGPaversym,VGPint=3)
          }else if(input$MVGP_Pole_Stat==4 && exists(MVGP_temp_stat_name, envir = MVGP_temp)==TRUE){
            plot_selected_VGP(s=s,lat0 = centerLat,lon0 = centerLong)
            plot_VGP_S(VGP = MVGP_temp[[MVGP_temp_stat_name]][[3]][,1:2],
                       lat = centerLat,
                       lon = centerLong,
                       auto_cent = F,
                       on_plot = TRUE, col = rgb(0,0,1,0.15),col_sym_out = rgb(0,0,1,0.15))
            PmagDiR::plot_PA95(lon = MVGP_temp[[MVGP_temp_stat_name]][[2]][1,2],
                               lat = MVGP_temp[[MVGP_temp_stat_name]][[2]][1,3],
                               A = 0,
                               lat0 = centerLat,
                               lon0 = centerLong,
                               col_f = col_f,symbol = MVGPaversym,size = 1.5,on_plot = T)
          }
          if(input$MVGP_Pole_Stat==4){
            if(exists(MVGP_temp_stat_name,envir = MVGP_temp)==FALSE){assign(x = MVGP_temp_stat_name,value = MVGP_boot_on_VGP,envir = MVGP_temp)}
            #extract result and assign name for export table
            MVGP_ALLVGPS_boot <- MVGP_temp[[MVGP_temp_stat_name]][[2]]
            rownames(MVGP_ALLVGPS_boot) <- input$fileN_MVGP

            #add name if requested
            if(input$MVGP_names_YN==2){
              N_MVGP <- input$fileN_MVGP
              #add name
              text(x = c2x(MVGP_temp[[MVGP_temp_stat_name]][[2]][1,2],MVGP_temp[[MVGP_temp_stat_name]][[2]][1,3]),
                   y = c2y(MVGP_temp[[MVGP_temp_stat_name]][[2]][1,2],MVGP_temp[[MVGP_temp_stat_name]][[2]][1,3]),
                   pos=3,
                   substitute(paste(bold(N_MVGP))), cex= 1.5)
            }
          }
        }
        #populate and plot CURRENT VGP pole statistic
        output$MVGP_ALLVGPS_stat <- renderTable({
          if(input$MVGP_Pole_Stat==1) {MVGP_ALLVGPS_stat <- NULL}
          if(input$MVGP_Pole_Stat==2) {MVGP_ALLVGPS_stat <- MVGP_ALL_Fisher}
          if(input$MVGP_Pole_Stat==3) {MVGP_ALLVGPS_stat <- MVGP_ALL_Fisher}
          if(input$MVGP_Pole_Stat==4) {MVGP_ALLVGPS_stat <- MVGP_ALLVGPS_boot}
          MVGP_ALLVGPS_stat
        }, digits = 1,align = "l", rownames = T)
      }

      #Add external Poles from list
      if(length(ext)){
        for(i in ext){
          PmagDiR::plot_PA95(lon = Added_poles$list[i,4],
                             lat = Added_poles$list[i,5],
                             A = Added_poles$list[i,6],
                             col_A = rgb(0,1,1,0.15),
                             lon0 = centerLong,
                             lat0 = centerLat,
                             col_f = Added_poles$list[i,2],
                             symbol = Added_poles$list[i,3],
                             on_plot = TRUE
          )
          if(input$addextreanames==2){
            polename <- Added_poles$list[i,1]
            #add name
            text(x = c2x(Added_poles$list[i,4],Added_poles$list[i,5]),
                 y = c2y(Added_poles$list[i,4],Added_poles$list[i,5]),
                 pos=3,
                 substitute(paste(bold(polename))), cex= 1.2)
          }
        }
        #Add Fisher of internal+manually added poles, if they are selected
        if(input$extrapolesfisher==2){
          #select color of average
          if(input$extrameansymbol==1) extrameansymbol <- "c"
          if(input$extrameansymbol==2) extrameansymbol <- "s"
          if(input$extrameansymbol==3) extrameansymbol <- "d"
          if(input$extrameansymbol==4) extrameansymbol <- "t"

          #select color of average
          if(input$extrameancolor==1) extrameancolor <- "black"
          if(input$extrameancolor==2) extrameancolor <- "blue"
          if(input$extrameancolor==3) extrameancolor <- "green"
          if(input$extrameancolor==4) extrameancolor <- "pink"
          if(input$extrameancolor==5) extrameancolor <- "purple"
          if(input$extrameancolor==6) extrameancolor <- "brown"
          if(input$extrameancolor==7) extrameancolor <- "red"
          if(input$extrameancolor==8) extrameancolor <- "yellow"
          if(input$extrameancolor==9) extrameancolor <- "cyan"
          if(input$extrameancolor==10) extrameancolor <- "gray"
          if(input$extrameancolor==11) extrameancolor <- "white"

          #add color and symbol to reactive file for saving in table
          extpole$extrameansymbol <- extrameansymbol
          extpole$extrameancolor <- extrameancolor

          #the extfisher reactive function is below!!!
          Ext_P_fisher <- extfisher()
          PmagDiR::plot_PA95(lon = Ext_P_fisher[1,1],
                             lat = Ext_P_fisher[1,2],
                             A = Ext_P_fisher[1,3],
                             lon0 = centerLong,
                             lat0 = centerLat,
                             col_A = rgb(1,0,0,0.2),
                             symbol = extrameansymbol,
                             col_f = extrameancolor,
                             size = 1.2,
                             on_plot = T)
          if(input$addextreanames==2){
            averpolename <- input$extreameanname
            #add name
            text(x = c2x(Ext_P_fisher[1,1],Ext_P_fisher[1,2]),
                 y = c2y(Ext_P_fisher[1,1],Ext_P_fisher[1,2]),
                 pos=3,
                 substitute(paste(bold(averpolename))), cex= 1.5)
          }
        }
      }
      #record plot
      MVGP_plot <- recordPlot()

      #Export graphic
      output$MVGP_G <- downloadHandler(
        filename = function() {
          paste(input$fileN_MVGP,"_Poles_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 11,height = 11)
          replayPlot(MVGP_plot)
          dev.off()
        }
      )
      output$MVGP_G2 <- downloadHandler(
        filename = function() {
          paste(input$fileN_MVGP,"_Poles_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 11,height = 11)
          replayPlot(MVGP_plot)
          dev.off()
        }
      )
      output$MVGP_G3 <- downloadHandler(
        filename = function() {
          paste(input$fileN_MVGP,"_Poles_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 11,height = 11)
          replayPlot(MVGP_plot)
          dev.off()
        }
      )


      #export Poles list and stat
      output$MVGP_list <- downloadHandler(
        filename = function() {
          paste(input$fileN_MVGP,"_Pole_list_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(MVGP_list$vgps,file,row.names = F)
        }
      )

      #export all VGPs stat
      #export locality level stat
      output$VGP_ALL_stat <- downloadHandler(
        filename = function() {
          paste(input$fileN_MVGP,"_POLE_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          if(input$MVGP_Pole_Stat==1){write.csv(NULL,file)}
          else if(input$MVGP_Pole_Stat==2){
            write.csv(round(MVGP_ALL_Fisher, digits = 2),file)
          }
          else if(input$MVGP_Pole_Stat==3){
            write.csv(round(MVGP_Fish_on_VGP, digits = 2),file)
          }
          else if(input$MVGP_Pole_Stat==4){
            write.csv(round(MVGP_ALLVGPS_boot, digits = 2),file)
          }
        }
      )
    }


    ############################################################

    ### Simulated VGP Part
    #function that generates vgp with Fisher pars
    SVGP_generator <- function(N,k,lon,lat,k_tol){
      #sub-function generating random long lat
      fisherDiR <- function(k){
        r2d <- function(x) {x*(180/pi)}
        L <- exp(-2*k)
        a <- runif(1)*(1-L)+L
        f <- sqrt(-log(a)/(2*k))
        latitude <- 90-r2d(2*asin(f))
        longitude <- r2d(2*pi*runif(1))
        return(c(longitude, latitude))
      }
      #reiterate until k is within the tolerance
      repeat{
        result <- data.frame(matrix(ncol = 2,nrow = 0))
        colnames(result) <- c("lon", "lat")
        for(i in 1:N){
          decinc_temp <- data.frame(t(fisherDiR(k)))
          colnames(decinc_temp) <- c("lon", "lat")
          result <- rbind(result,decinc_temp)
        }
        AverLongLat <- PmagDiR::fisher(result)
        fixed_data <- PmagDiR::bed_DI(result,in_file = F,
                                      bed_az = (AverLongLat[1,1]+180)%%360,bed_plunge =90-AverLongLat[1,2])
        final_VGP <- PmagDiR::bed_DI(fixed_data,in_file = F,
                                     bed_az = lon,bed_plunge =90-lat)
        colnames(final_VGP) <- c("Long","Lat")
        stat <- PmagDiR::fisher(final_VGP)
        k_test <- stat[1,6]
        #check for tolerance
        if(k_tol==0) break
        if(abs(k_test-k)<=k_tol) break
      }
      return(final_VGP)
    }

    #creates reactive files
    SVGP <- reactiveValues(vgps = NULL)
    SVGP <- reactiveValues(vgps_old = NULL)
    SVGPS <- reactiveValues(stat = NULL)

    #populate reactive file
    observeEvent(input$SVGPgo,{
      SVGP$vgps <- SVGP_generator(N = input$SVGPN,k = input$SVGPk,lon = input$SVGPlon,lat = input$SVGPlat,k_tol=input$k_tol)
    })

    #function that generate plot
    plot_SVGP <- function(){
      req(SVGP$vgps)
      VGP <- SVGP$vgps

      #choose color
      if(input$SVGPcolor==1) SVGPcolor <- "black"
      if(input$SVGPcolor==2) SVGPcolor <- "blue"
      if(input$SVGPcolor==3) SVGPcolor <- "green"
      if(input$SVGPcolor==4) SVGPcolor <- "pink"
      if(input$SVGPcolor==5) SVGPcolor <- "purple"
      if(input$SVGPcolor==6) SVGPcolor <- "brown"
      if(input$SVGPcolor==7) SVGPcolor <- "red"
      if(input$SVGPcolor==8) SVGPcolor <- "yellow"
      if(input$SVGPcolor==9) SVGPcolor <- "cyan"
      if(input$SVGPcolor==10) SVGPcolor <- "gray"
      if(input$SVGPcolor==11) SVGPcolor <- "white"
      assign("SVGP_color",SVGPcolor, envir = MVGP_temp)

      #choose symbol
      if(input$SVGPsymbol==1) SVGPsymbol <- "c"
      if(input$SVGPsymbol==2) SVGPsymbol <- "s"
      if(input$SVGPsymbol==3) SVGPsymbol <- "d"
      if(input$SVGPsymbol==4) SVGPsymbol <- "t"
      assign("SVGP_symbol",SVGPsymbol, envir = MVGP_temp)

      SVGPS$stats <- plot_VGP_S(VGP = VGP,lat = input$SVGP_clat,long = input$SVGP_clong,
                                coast = ifelse(input$SVGP_coast==1,TRUE,FALSE),
                                auto_cent = ifelse(input$SVGP_center==1,TRUE,FALSE),
                                col = SVGPcolor,symbol = SVGPsymbol,nb = input$SVGPnb,
                                A95 = ifelse(input$SVGP_stat==2,TRUE,FALSE),
                                B95 = ifelse(input$SVGP_stat==3,TRUE,FALSE),
                                # B95 = bootTF,
                                VGPint = 4)


      output$SVGPstat <- renderTable({
        if(input$SVGP_stat==1) {SVGPpolestat <- NULL}
        if(input$SVGP_stat==2) {
          SVGPpolestat <- SVGPS$stats[[1]]
          rownames(SVGPpolestat) <- input$SVGPname
        }
        if(input$SVGP_stat==3) {
          SVGPpolestat <- SVGPS$stats[[2]]
          rownames(SVGPpolestat) <- input$SVGPname
        }
        SVGPpolestat
      }, digits = 1,align = "l", rownames = T)
    }

    #send current VGP to list and reactive VGP_saved_list
    observeEvent(input$saveSVGP,{
      req(SVGP$vgps)
      VGP <- SVGP$vgps
      colnames(VGP) <- c("Long","Lat")
      VGP$color <- MVGP_temp[["SVGP_color"]]
      VGP$symbol <- MVGP_temp[["SVGP_symbol"]]
      VGPsaved <- list()
      VGPsaved[[1]] <- VGP
      #if statistics does not exists, it makes it
      if(length(SVGPS$stats)<=1){
        VGPsaved[[2]] <- plot_VGP_S(VGP = VGP,A95 = T,B95 = T,nb = input$SVGPnb,VGPint = 4,plot=F)
      }else{VGPsaved[[2]] <- SVGPS$stats}
      VGP_saved$list[[input$SVGPname]] <- VGPsaved
      MVGP_list$vgps <- sites_list()
    })

    #preliminary table for loaded VGP
    output$Ext_VGP_list2 <- renderTable({
      if(length(VGP_saved$list)==0){
        MVGP_list$vgps <- NULL
      }else{
        MVGP_list$vgps <- sites_list()
      }
      MVGP_list$vgps
    }, digits = 1,align = "l", rownames = T)

    #send to ui
    output$SVGP_plot <- renderPlot({
      plot_SVGP()
      #record plot
      SVGP_plot <- recordPlot()

      #Export graphic
      output$SVGP_G <- downloadHandler(
        filename = function() {
          paste(input$SVGPname,"_Siumlated_VGP_N_",input$SVGPN,"_K_",input$SVGPk,"_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 11,height = 11)
          replayPlot(SVGP_plot)
          dev.off()
        }
      )
      #export stat
      output$SVGP_S <- downloadHandler(
        filename = function() {
          paste(input$SVGPname,"_Pole_Siumlated_VGP_N_",input$SVGPN,"_K_",input$SVGPk,"_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          if(input$SVGP_stat==1){write.csv(NULL,file)}
          else if(input$SVGP_stat==2){
            SVGP_pole <- SVGPS$stats[[1]]
            rownames(SVGP_pole) <- paste(input$SVGPname,"_fisher",sep="")
            write.csv(round(SVGP_pole, digits = 2),file)
          }
          else if(input$SVGP_stat==3){
            SVGP_pole <- SVGPS$stats[[2]]
            rownames(SVGP_pole) <- paste(input$SVGPname,"_bootstrap",sep="")
            write.csv(round(SVGP_pole, digits = 2),file)
          }
        }
      )
      output$SVGP_list <- downloadHandler(
        filename = function() {
          paste(input$SVGPname,"_List_Siumlated_VGP_N_",input$SVGPN,"_K_",input$SVGPk,"_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(round(SVGP$vgps,digits = 2),file = file,row.names = F)
        }
      )
    }, width = 700,height = 700)


    ####ADD EXTERNAL POLES PART
    #create reactive tab file
    Added_poles <- reactiveValues(list=NULL)

    #creates reactive value for checking if file is uploaded
    values <- reactiveValues(mvgp = NULL)
    #check for uploaded file
    observeEvent(input$vgpfile,{values$mvgp <- "uploaded"})
    #reset upload if requested
    observeEvent(input$resetMVGP,{values$mvgp <- "reset"})
    #read file if present, reset if requested
    VGPfile <- reactive({
      if (is.null(values$mvgp)) {
        return(NULL)
      } else if (values$mvgp == 'uploaded') {
        read.csv(file = input$vgpfile$datapath)
      } else if (values$mvgp == 'reset') {
        return(NULL)
      }
    })

    #creates reactive value for checking if file is uploaded
    extpole <- reactiveValues(listfile = NULL)
    #check for uploaded file
    observeEvent(input$extrapolelist,{extpole$listfile <- "uploaded"})
    #reset upload if requested
    observeEvent(input$delextrapolelist,{extpole$listfile <- "reset"})
    #read file if present, reset if requested
    Extra_poles_list <- reactive({
      if(is.null(extpole$listfile)){
        return(NULL)
      }else if(extpole$listfile == "uploaded"){
        temp <- read.csv(file = input$extrapolelist$datapath)
        colnames(temp) <- c("Loc.","Col","Sym","Lon","Lat","A95")
        return(temp)
      }else if(extpole$listfile == "reset"){
        return(NULL)
      }
    })

    #add manual pole to list if exists
    observeEvent(input$addextrapole,{
      empty_list <- data.frame(matrix(ncol=6,nrow = 0))
      colnames(empty_list) <- c("Loc.","Col","Sym","Lon","Lat","A95")
      Added_poles$list <- empty_list
      #choose color of average
      if(input$extrapolecolor==1) EPcol <- "black"
      if(input$extrapolecolor==2) EPcol <- "blue"
      if(input$extrapolecolor==3) EPcol <- "green"
      if(input$extrapolecolor==4) EPcol <- "pink"
      if(input$extrapolecolor==5) EPcol <- "purple"
      if(input$extrapolecolor==6) EPcol <- "brown"
      if(input$extrapolecolor==7) EPcol <- "red"
      if(input$extrapolecolor==8) EPcol <- "yellow"
      if(input$extrapolecolor==9) EPcol <- "cyan"
      if(input$extrapolecolor==10) EPcol <- "gray"
      if(input$extrapolecolor==11) EPcol <- "white"

      #select symbol of average
      if(input$extrapolesimbol==1) EPsym <- "c"
      if(input$extrapolesimbol==2) EPsym <- "s"
      if(input$extrapolesimbol==3) EPsym <- "d"
      if(input$extrapolesimbol==4) EPsym <- "t"

      temp <- data.frame(matrix(ncol=6,nrow = 0))
      colnames(temp) <- c("Loc.","Col","Sym","Lon","Lat","A95")
      temp[1,1] <- input$extrapolename
      temp[1,2] <- EPcol
      temp[1,3] <- EPsym
      temp[1,4] <- input$extrapolelong
      temp[1,5] <- input$extrapolelat
      temp[1,6] <- input$extrapoleA95

      if(is.na(input$extrapolelong)==F && is.na(input$extrapolelat)==F && is.na(input$extrapoleA95)==F) {
        Added_poles$list <- rbind(Added_poles$list,temp)
      }
      if(!is.null(extpole$listfile)){Added_poles$list <- rbind(Added_poles$list,Extra_poles_list())}

      #eliminates duplicates
      Added_poles$list <- dplyr::distinct(Added_poles$list)
      if(nrow(Added_poles$list)==0){Added_poles$list <- NULL}
    })

    #delete elements from list
    observeEvent(input$deleteextrapole,{
      d <- input$EP_list_rows_selected
      if(length(d)){Added_poles$list <- Added_poles$list[-d,]}
      if(nrow(Added_poles$list)==0){Added_poles$list <- NULL}
    })

    #send to ui the interactive extra pole list
    output$EP_list <- DT::renderDataTable(Added_poles$list, server = F)

    #Calculate Fisher of selected poles, SENT TO UI within "all_poles_plotter"

    extfisher <- reactive({
      if(input$extrapolesfisher==2){
        extlist <- data.frame(matrix(ncol=2,nrow=0))
        colnames(extlist) <- c("Lon","Lat")

        s <- input$MVGPlist_rows_selected
        ext <- input$EP_list_rows_selected
        if(length(s)){
          for(i in s){
            lonlat_temp <- MVGP_list$vgps[i,5:6]
            lonlat_temp$Lon <- as.numeric(lonlat_temp$Lon)
            lonlat_temp$Lat <- as.numeric(lonlat_temp$Lat)
            extlist <- rbind(extlist,lonlat_temp)
          }
        }
        if(length(ext)){
          for(e in ext){
            lonlat_temp <- Added_poles$list[e,4:5]
            extlist <- rbind(extlist,lonlat_temp)
          }
        }
        if(length(s)!=0 || length(ext)!=0){
          Fisher <- PmagDiR::fisher(extlist)
          Fisher <- Fisher[,-5]
          colnames(Fisher) <- c("Long","Lat","A95","N","K")
          rownames(Fisher) <- input$extreameanname
        }else{Fisher <- NULL}
        extpole$Fisher <- Fisher
        return(Fisher)
      }
    })
    #send fisher to UI
    output$extpolesfisher <- renderTable({
      if(input$extrapolesfisher==2){
        extpole$Fisher
      }else{NULL}
    },digits = 1,align = "l", rownames = T)

    #send fisher to download
    output$Ext_pole_fisher_S <- downloadHandler(
      filename = function(){
        paste(input$fileN_MVGP,"_All_Poles_Fisher_", Sys.Date(), ".csv", sep="")
      },
      content = function(file){
        Ext_poles_S <- extfisher()
        write.csv(round(Ext_poles_S,digits = 2),file = file,row.names = T)
      }
    )

    #Add external pole list Fisher to external pole list
    observeEvent(input$add_F_2_ext_list,{
      if(is.null(extpole$Fisher)==F){
        temp <- data.frame(matrix(ncol=6,nrow = 1))
        colnames(temp) <- c("Loc.","Col","Sym","Lon","Lat","A95")
        temp[1,1] <- input$extreameanname
        temp[1,2] <- extpole$extrameancolor
        temp[1,3] <- extpole$extrameansymbol
        temp[1,4:6] <- round(extpole$Fisher[1,1:3],digits = 2)
        Added_poles$list <- rbind(Added_poles$list,temp)
      }
    })


    #############SEND PLOTS TO VGP ANALYSIS FIGURES

    #send figure to ui of Multi VGP
    output$MVGP_plot <- renderPlot({all_poles_plotter()},width = 700, height = 700)

    #send figure to ui of Multi VGP2 equal as above
    output$MVGP_plot2 <- renderPlot({all_poles_plotter()},width = 700, height = 700)

    #send figure to ui of Multi VGP3 equal as above
    output$MVGP_plot3 <- renderPlot({all_poles_plotter()},width = 700, height = 700)


    ############ END OF VIRTUAL GEOMAGNETIC POLES MODULE

    ############ UNSTRAIN MODULE
    #read matrix file
    Str_M <- reactive({
      if(is.null(input$str_mtrx)){
        return("")
      }
      read.csv(file = input$str_mtrx$datapath)
    })

    #table with inverted AMS eigenvectors
    output$straindirs <- renderTable({
      #select type of matrix to invert
      if(input$str_m_type==1){type <- "v"}
      else if(input$str_m_type==2){type <- "m"}
      else if(input$str_m_type==3){type <- "m3x3"}

      #invert matrix
      if(is.null(input$str_mtrx)==FALSE){
        inv_result <- AMS_inv(Str_M(),type = type,prnt=F, Shiny = TRUE)
        assign("inv_result",inv_result,envir = .GlobalEnv)
        tab <- inv_result[[2]]
        tab
      }
    },rownames = T, caption="Eigenvectors (dec and inc of V1,V2,V3), Lineation (L_inv), and Foliation (F_inv) of the inverse strain matrix")

    #unstrain function modified for Web-PmagDiR
    unstrain <- eventReactive(input$unstr_GO,{
      if(exists("inv_result", envir = .GlobalEnv)==FALSE){stop("LOAD STRAIN PARAMETER BEFORE EXECUTE")}
      S_vec <- inv_result[[1]]
      #unstrain function adapted for Shiny
      unstr_DI_S <- function(DIAP,S_vec,Lin,Fol,n=1,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE){
        #degree to radians and vice versa
        d2r <- function(x) {x*(pi/180)}
        r2d <- function(x) {x*(180/pi)}

        #load directions data
        data <- DIAP
        data <- na.omit(data)
        colnames(data) <- c("dec", "inc","baz","binc")
        dirs <- data[,1:2]
        bed <- data[,3:4]

        #directions in Cartesian coordinates
        dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
        dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
        dirs$z <- sin(d2r(dirs$inc))

        #bedding in Cartesian coordinates
        bed$x <- cos(d2r(bed$baz))*cos(d2r(bed$binc))
        bed$y <- sin(d2r(bed$baz))*cos(d2r(bed$binc))
        bed$z <- sin(d2r(bed$binc))

        inc_e_dec <- as.data.frame(matrix(ncol=3,nrow = 1))
        BdecInc <- bed_DI(data[,1:4])
        inc_e_dec[1,1:3] <- as.data.frame(inc_E_finder(BdecInc))
        inc_e_dec[1,3] <- abs(inc_e_dec[1,3])
        colnames(inc_e_dec) <- c("V1inc","E","DV1V2")

        #set parameters of deforming matrix
        Lincr <- (Lin-1)/n
        Fincr <- (Fol-1)/n
        L <- 1
        F <- 1
        ninc <- 0
        repeat{
          ninc <- ninc+1
          #lineation
          L <- L+Lincr
          #foliation
          F <- F+Fincr
          #anisotropy degree
          P <- F*L
          #eigenvalues
          K1 <- (3*L)/(L+(1/F)+1)
          K2 <- K1/L
          K3 <- K1/P

          #matrix of new eigenvalue
          M <- c(K1,0,0,0,K2,0,0,0,K3)
          M <- matrix(M,nrow=3,byrow=T)

          #combines given eigenvalues with Strain directions
          S <- S_vec%*%M%*%inv(S_vec)
          S_e <- eigen(S,symmetric = T)
          S_vec <- S_e$vectors
          S_val <- S_e$values
          new_DI <- data.frame(matrix(ncol=4,nrow=0))
          for (i in 1:nrow(data)){
            #unstrain dirs
            dircart <- t(as.matrix(dirs[i,3:5]))
            unstr <- S%*%dircart
            NewInc <- r2d(asin(unstr[3,1]/(sqrt((unstr[1,1]^2)+(unstr[2,1]^2)+(unstr[3,1]^2)))))
            NewDec <- r2d(atan2(unstr[2,1],unstr[1,1]))
            NewDec <- NewDec%%360
            #unstrain bedding
            bedcart <- t(as.matrix(bed[1,3:5]))
            bunstr <- S%*%bedcart
            Newbinc <- r2d(asin(bunstr[3,1]/(sqrt((bunstr[1,1]^2)+(bunstr[2,1]^2)+(bunstr[3,1]^2)))))
            Newbaz <- r2d(atan2(bunstr[2,1],bunstr[1,1]))
            Newbaz <- Newbaz%%360

            new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc),
                              as.data.frame(Newbaz),as.data.frame(Newbinc))
            new_DI <- rbind(new_DI,new_DI_p)
          }
          colnames(new_DI) <- c("dec","inc","baz","binc")
          NBdecInc <- bed_DI(new_DI)
          inc_e_dec_p <- inc_E_finder(NBdecInc)
          inc_e_dec <- rbind(inc_e_dec,inc_e_dec_p)
          inc_e_dec[,3] <- abs(inc_e_dec[,3])

          #generate tk03 curve also for later plot
          tkx <- 0:90
          tky <- tk03(tkx)

          #function that update the progress bar of shiny
          updateProgressBar(
            id="nincrements",
            value=ninc,total=input$increm,
          )

          if (cross==TRUE){
            #break loop if crosses tk03 line
            if(any(inc_e_dec$E>tk03(inc_e_dec$V1inc)) && any(inc_e_dec$E<tk03(inc_e_dec$V1inc))){
              curve1 <- inc_e_dec[,1:2]
              curve1 <- curve1[order(curve1$E),]
              curve2 <- cbind(as.data.frame(tkx),as.data.frame(tky))
              crxy <- curve_cross(curve1,curve2)
              break
            }
          }
          #break loops when maximum Edec reached
          if(EdMAX==TRUE && length(inc_e_dec[,3])>3){
            l <- length(inc_e_dec[,3])
            if(inc_e_dec[l,3]<inc_e_dec[l-1,3] && inc_e_dec[l-1,3]>inc_e_dec[l-2,3]){
              inc_e_dec <- inc_e_dec[-l,]
              break
            }
          }
          #break loops when minimu Edec reached
          if(EdMIN==TRUE && length(inc_e_dec[,3])>3){
            l <- length(inc_e_dec[,3])
            if(inc_e_dec[l,3]>inc_e_dec[l-1,3] && inc_e_dec[l-1,3]<inc_e_dec[l-2,3]){
              inc_e_dec <- inc_e_dec[-l,]
              break
            }
          }
          #break loop if reaches defined F and L
          if(round(L,digits = 4) == round(Lin,digits = 4) && round(F,digits = 4) == round(Fol,digits = 4)) {
            break
          }
          finalpars <- paste("Final L: ",round(L,digits = 5), ", Final F:",round(F,digits = 5))
          assign("finalpars",finalpars,envir = .GlobalEnv)
        }
        #set figure
        par(fig=c(0,1,0,1), new= FALSE)
        y_up <- ifelse(inc_e_dec[1,2]>3.5, ceiling(inc_e_dec[1,2]),3.5)

        plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
             xlab="Inclination (°)", ylab="Elongation", cex.lab=1.5)

        #text for figure
        N <- as.character(length(data[,1]))
        Inc <- round(inc_e_dec[1,1],digits = 1)
        Ecut <- round(inc_e_dec[1,2],digits = 2)
        V2 <- round(inc_e_dec[1,3],digits = 1)
        text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
        text(x=0, y=3.2,pos=4,text, cex= 1.2)

        #plot tk03.GAD model E-I
        points(x=tkx, y= tky, type= "l", col="blue", lwd=3)
        #plot unstrain curve
        points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="red", lwd=3)
        points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
               pch=21, col="black", bg="blue", cex=1.4)
        #plot cross with EI line if exists
        if(exists("crxy")==TRUE){
          points(x=crxy[1,1],y=crxy[1,2],
                 pch=21, col="black", bg="red", cex=1.4)
        }else{
          points(x=inc_e_dec[length(inc_e_dec[,1]),1],
                 y=inc_e_dec[length(inc_e_dec[,1]),2],
                 pch=21, col="black", bg="red", cex=1.4)
        }
        #unstrained parameters for text2 of figure
        inc_nstr <- ifelse(exists("crxy")==TRUE,
                           round(crxy[1,1], digits=1),
                           round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1))

        E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
        Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

        #put text unstrained in figure
        text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
        text(x=20, y=3.2, pos=4, text2, cex= 1.2)

        #plot Edec during unstrain
        par(fig=c(0.6,1,0.56,0.99), new=TRUE)
        plot(x=inc_e_dec$V1inc, y=inc_e_dec$DV1V2, type="l", tck=-0.05,
             frame.plot=FALSE,cex.axis=0.8, mgp=c(1.5,0.5,0.1),
             xlab="",ylab="",yaxt="n",
             ylim=c(floor(min(inc_e_dec$DV1V2)),
                    ceiling(max(inc_e_dec$DV1V2))),
             xlim=c(floor(min(inc_e_dec$V1inc)),
                    ceiling(max(inc_e_dec$V1inc))))
        axis(2, cex.axis=0.8, las=2)
        title(xlab = "Inc (°)", line=1.6)
        title(ylab = "Edec (°)", line=2.5)

        points(x=inc_e_dec[1,1],y=inc_e_dec[1,3],
               pch=21, col="black", bg="blue", cex=1.4)
        points(x=inc_e_dec[length(inc_e_dec[,1]),1],
               y=inc_e_dec[length(inc_e_dec[,1]),3],
               pch=21, col="black", bg="red", cex=1.4)

        #sets results file
        unstr_matrix <- S
        unstr_results <- list(data,BdecInc,NBdecInc, new_DI,inc_e_dec,unstr_matrix)
        names(unstr_results) <- c("Original dataset","original TC directions","unstrained TC directions",
                                  "unstrained directions and bedding",
                                  "inc, E, declination triplets","Unstrain matrix")
        return(unstr_results)
      }
      #performe routine
      Unstrain_Dirs <- function(){
        DIAP <- fix_DI(input_file(),coord = 1)
        if(ncol(DIAP)<4) stop("File must contain declination, inclination and bedding in 4 columns")
        unstr_DI_S(DIAP = DIAP,S_vec = S_vec,Lin = input$lin,Fol = input$fol,n=input$increm,
                   cross = ifelse(input$brk==2, TRUE, FALSE),
                   EdMAX = ifelse(input$brk==3, TRUE, FALSE),
                   EdMIN = ifelse(input$brk==4, TRUE, FALSE))
      }
      unstr_result <- Unstrain_Dirs()
      assign("unstr_result",unstr_result,envir = .GlobalEnv)
      #paste text final pars
      output$finalpars <- renderText(finalpars)
      #plot corrected directions, always tilt corrected coords
      output$unstrDirs <- renderPlot({

        #fix unstrained dirs coordinate depending on input
        fix_DI_str <- function(DIAP){
          req(input$lat)
          req(input$long)
          Slat <- input$lat
          Slong <- input$long
          if(input$coord_str==1){DI <- DIAP}
          else if(input$coord_str==2){DI <- bed_DI(DIAP)}

          #apply cutoff
          if(input$coord_str==1){
            geo=TRUE
          }else{geo=FALSE}
          if(input$cutoff_str==2){DI <- PmagDiR::cut_DI(DI = DI,lat=Slat,long = Slong,geo = geo,Shiny = T)}
          else if(input$cutoff_str==3){DI <- PmagDiR::cut_DI(DI = DI,lat=Slat,long = Slong,inc_f = F,geo = geo,Shiny = T)}
          else if(input$cutoff_str==4){DI <- PmagDiR::cut_DI(DI = DI,VD=F,cutoff = input$VGP_fixed_str ,lat=Slat,long = Slong,geo = geo,Shiny = T)}
          else if(input$cutoff==5){DI <- PmagDiR::cut_DI(DI = DI,VD=F,cutoff = input$VGP_fixed_str ,lat=Slat,long = Slong, inc_f=F,geo = geo,Shiny = T)}
          else {DI <- DI}
          return(DI)
        }

        #equal area function
        plot_EA <- function(dirs){
          req(input$lat)
          req(input$long)
          Slat <- input$lat
          Slong <- input$long
          DI <- fix_DI_str(dirs)
          if(input$mode_str==1){DI <- DI}
          if(input$mode_str==2){DI <- common_DI(DI)}
          if(input$mode_str==3){DI <- common_DI(DI,down = F)}
          #next just flip negative to positive or vice versa when required
          if(input$mode_str==4){
            for(i in 1:nrow(DI)){
              if(DI[i,2]<0){
                DI[i,1] <- (DI[i,1]+180)%%360
                DI[i,2] <- abs(DI[i,2])
              }
            }
          }else if(input$mode_str==5){
            for(i in 1:nrow(DI)){
              if(DI[i,2]>=0){
                DI[i,1] <- (DI[i,1]+180)%%360
                DI[i,2] <- -(DI[i,2])
              }
            }
          }
          #define colors Down-pointing
          if(input$colD_str==1) colD <- "black"
          if(input$colD_str==2) colD <- "blue"
          if(input$colD_str==3) colD <- "red"
          if(input$colD_str==4) colD <- "dark green"

          #define color Up-pointing
          if(input$colU_str==1) colU <- "white"
          if(input$colU_str==2) colU <- "cyan"
          if(input$colU_str==3) colU <- "pink"
          if(input$colU_str==4) colU <- "light green"

          #define symbol
          if(input$sym_str==1) sym <- "c"
          if(input$sym_str==2) sym <- "s"
          if(input$sym_str==3) sym <- "d"
          if(input$sym_str==4) sym <- "t"

          plot_DI(DI,col_d = colD,col_u = colU, symbol = sym)
          if(input$fisher_str==2){
            unstr_DI_stat <- fisher_plot_S(DI)
          }else if(input$fisher_str==3){
            unstr_DI_stat <- ellips_plot_S(DI,lat = Slat,long = Slong)
          }else{unstr_DI_stat <- NULL}
          #swrite table
          output$unstr_DI_stat <- renderTable(unstr_DI_stat, rownames = T)
          output$unstrdirsS <- downloadHandler(
            filename = function() {
              paste(input$fileN_str,"_unstr_DI_", Sys.Date(), "_stat.csv", sep="")
            },
            content = function(file) {
              unstr_DI_stat <- round(unstr_DI_stat, digits = 2)
              write.csv(unstr_DI_stat, file)
            }
          )
        }
        #execute funtion
        plot_EA(unstr_result[[4]])
        #replot for download
        output$unstrdirsG <- downloadHandler(
          filename = function(){
            paste(input$fileN_str,"_unstr_DI_", Sys.Date(), ".pdf", sep="")
          },
          content = function(file){
            pdf(file, onefile = TRUE,width = 9,height = 9)
            plot_EA(unstr_result[[4]])
            dev.off()
          }
        )
        output$unstrdirsD <- downloadHandler(
          filename = function(){
            paste(input$fileN_str,"_unstr_DI_", Sys.Date(), ".csv", sep="")
          },
          content = function(file){
            DI <- fix_DI_str(unstr_result[[4]])
            write.csv(round(DI,digits = 2),file,row.names = F)
          }
        )


      },width = 700,height = 700)
      #save plot for export
      unstr_fig <- recordPlot()
      #export plot
      output$unstrG <- downloadHandler(
        filename = function() {
          paste(input$fileN_str,"_unstr_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 12,height = 11)
          replayPlot(unstr_fig)
          dev.off()
        }
      )


    })
    #unstrain bootstrap function for web-PmagDiR
    unstrain_b <- eventReactive(input$unstr_boot_GO, {

      #unstrain bootstrap function adapted for Shiny
      unstr_boot_S <- function(unstr_file,nb= 100,S_vec,Lin,Fol,ns=1,confidence=95,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE,hist=TRUE){

        #degree to radians and vice versa
        d2r <- function(x) {x*(pi/180)}
        r2d <- function(x) {x*(180/pi)}

        #load directions data
        dat <- unstr_file[[1]]
        dat <- na.omit(dat)
        colnames(dat) <- c("dec", "inc","baz","binc")

        #load inc_e_dec result
        inc_e_dec <- unstr_file[[5]]

        #plot frame
        par(fig=c(0,1,0,1), new= FALSE)
        y_up <- ifelse(inc_e_dec[1,2]>3.5, ceiling(inc_e_dec[1,2]),3.5)

        plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
             xlab="Inclination (°)", ylab="Elongation",cex.lab=1.5)

        #plot tk03.GAD model E-I
        x <- 0:90
        y <- tk03(x)
        points(x=x, y= y, type= "l", col="blue", lwd=3)

        #plot unstrain curve
        points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="red", lwd=3)
        points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
               pch=21, col="black", bg="blue", cex=1.2)
        points(x=inc_e_dec[length(inc_e_dec[,1]),1],y=inc_e_dec[length(inc_e_dec[,1]),2],
               pch=21, col="black", bg="red", cex=1.2)

        #text for figure
        N <- nrow(dat)
        Inc <- round(inc_e_dec[1,1],digits = 1)
        Ecut <- round(inc_e_dec[1,2],digits = 2)
        Edec <- round(inc_e_dec[1,3],digits = 1)
        text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",Edec)
        text(x=0, y=3.2,pos=4,text, cex= 1.2)

        #unstrained parameters for text2 of figure
        inc_nstr <- round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1)
        E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
        Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

        #put text unstrained in figure
        text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
        text(x=20, y=3.2, pos=4, text2, cex=1.2)

        #create empty files
        init_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
        final_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
        colnames(init_E_I) <- c("Inc","E","E_dec")
        colnames(final_E_I)<- c("Inc","E","E_dec")

        #start bootstrapping
        n <- 0
        repeat{
          n <- n+1
          Seq_I_E_b <- as.data.frame(matrix(ncol=3,nrow=0))
          E_declin <- as.data.frame(matrix(ncol=1,nrow=0))
          data <- boots_DI(dat)
          dirs <- data[,1:2]
          bed <- data[,3:4]

          #directions in Cartesian coordinates
          dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
          dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
          dirs$z <- sin(d2r(dirs$inc))

          #bedding in Cartesian coordinates
          bed$x <- cos(d2r(bed$baz))*cos(d2r(bed$binc))
          bed$y <- sin(d2r(bed$baz))*cos(d2r(bed$binc))
          bed$z <- sin(d2r(bed$binc))

          #set parameters of deforming matrix
          Lincr <- (Lin-1)/ns
          Fincr <- (Fol-1)/ns
          L <- 1
          F <- 1

          #gradual unstrain of the n pseudosample
          repeat{
            #anisotropy degree
            P <- F*L
            #eigenvalues
            K1 <- (3*L)/(L+(1/F)+1)
            K2 <- K1/L
            K3 <- K1/P

            #matrix of new eigenvalue
            M <- c(K1,0,0,0,K2,0,0,0,K3)
            M <- matrix(M,nrow=3,byrow=T)

            #combines given eigenvalues with Strain directions
            S <- S_vec%*%M%*%inv(S_vec)

            new_DI <- data.frame(matrix(ncol=4,nrow=0))
            for (i in 1:nrow(data)){
              #unstrain dirs
              dircart <- t(as.matrix(dirs[i,3:5]))
              unstr <- S%*%dircart
              NewInc <- r2d(asin(unstr[3,1]/(sqrt((unstr[1,1]^2)+(unstr[2,1]^2)+(unstr[3,1]^2)))))
              NewDec <- r2d(atan2(unstr[2,1],unstr[1,1]))
              NewDec <- NewDec%%360

              #unstrain bedding
              bedcart <- t(as.matrix(bed[i,3:5]))
              bunstr <- S%*%bedcart
              Newbinc <- r2d(asin(bunstr[3,1]/(sqrt((bunstr[1,1]^2)+(bunstr[2,1]^2)+(bunstr[3,1]^2)))))
              Newbaz <- r2d(atan2(bunstr[2,1],bunstr[1,1]))
              Newbaz <- Newbaz%%360

              new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc),
                                as.data.frame(Newbaz),as.data.frame(Newbinc))
              new_DI <- rbind(new_DI,new_DI_p)
            }
            colnames(new_DI) <- c("dec","inc","baz","binc")
            NBdecInc <- bed_DI(new_DI)
            inc_e_dec_p <- inc_E_finder(NBdecInc)
            Seq_I_E_b <- rbind(Seq_I_E_b,inc_e_dec_p)
            #take only absolute value of declination for the breaking conditions
            E_declin <- rbind(E_declin,abs(inc_e_dec_p[1,3]))

            #if cross== TRUE breaks loop when tk03.GAD is crossed
            if(cross==TRUE)  {if(any(Seq_I_E_b$E>tk03(Seq_I_E_b$V1inc))
                                 && any(Seq_I_E_b$E<tk03(Seq_I_E_b$V1inc))) break}

            #if EdMAX==TRUE breaks loop when max of Edec is reached
            if(EdMAX==TRUE && length(E_declin[,1])>3){
              l <- length(E_declin[,1])
              if(E_declin[l,1]<E_declin[l-1,1] && E_declin[l-1,1]>E_declin[l-2,1]){
                E_declin <- E_declin[-l,]
                break
              }
            }
            #break loops when minimu Edec reached
            if(EdMIN==TRUE && length(E_declin[,1])>3){
              l <- length(E_declin[,1])
              if(E_declin[l,1]>E_declin[l-1,1] && E_declin[l-1,1]<E_declin[l-2,1]){
                E_declin <- E_declin[-l,]
                break
              }
            }

            #break loop if reaches defined F and L
            if(round(L,digits = 4) == round(Lin,digits = 4) &&
               round(F,digits = 4) == round(Fol,digits = 4)) break
            #lineation
            L <- L+Lincr
            #foliation
            F <- F+Fincr
          }

          #if cross==TRUE select only curves that cross tk03.GAD
          if(cross==TRUE){
            if(any(Seq_I_E_b$E>tk03(Seq_I_E_b$V1inc)) && any(Seq_I_E_b$E<tk03(Seq_I_E_b$V1inc))==TRUE){
              points(x=Seq_I_E_b$V1inc, y= Seq_I_E_b$E,
                     type= "l", col=rgb(1, 0, 0, 0.15), lwd=1)
              i_E_I <- Seq_I_E_b[1,]
              f_E_I <- Seq_I_E_b[length(Seq_I_E_b[,1]),]
            }
            #if cross== FALSE, use any other condition to fill files
          }else{
            points(x=Seq_I_E_b$V1inc, y= Seq_I_E_b$E,
                   type= "l", col=rgb(0, 0.5, 1, 0.15), lwd=1)
            i_E_I <- Seq_I_E_b[1,]
            f_E_I <- Seq_I_E_b[length(Seq_I_E_b[,1]),]
          }
          colnames(i_E_I) <- c("Inc","E","E_dec")
          colnames(f_E_I)<- c("Inc","E","E_dec")
          init_E_I <- rbind(init_E_I,i_E_I)
          final_E_I <- rbind(final_E_I,f_E_I)
          init_E_I <- na.omit(init_E_I)
          final_E_I <- na.omit(final_E_I)

          #function that update the progress bar of shiny
          updateProgressBar(
            id="str_boots",
            value=nrow(final_E_I),total=input$strboot
          )


          if(nrow(final_E_I)==nb) {break}

        } #end of bootstrap

        #backup_files
        final_E_Ibk <- final_E_I
        init_E_Ibk <- init_E_I

        #Plot Bootstrapped data
        if(cross==FALSE){
          points(x=final_E_I$Inc,
                 y=final_E_I$E,
                 pch=16,
                 col=rgb(1,0,0,0.2),
                 cex=0.8)
        }

        #cut bootstrapped results for 95% (unless different) confidence
        conf <- confidence/100
        num <- round((nb*(1-conf))/2,digits=0)
        Lconf <- num
        Uconf <- nb-num
        final_E_I <- final_E_I[order(final_E_I$Inc),]
        final_E_I <- final_E_I[Lconf:Uconf,]

        points(x=x, y= y, type= "l", col="blue", lwd=3)

        #plot unstrain curve
        points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="yellow", lwd=2.2)
        points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
               pch=21, col="black", bg="blue", cex=1.5)
        #plot cross with EI line
        points(x=inc_e_dec[length(inc_e_dec[,1]),1],y=inc_e_dec[length(inc_e_dec[,1]),2],
               pch=21, col="black", bg=rgb(1,0.4,0.4,1), cex=1.5)

        #draw two lines for 95% confidence margin
        arrows(x0=final_E_I[1,1],x1=final_E_I[1,1],
               y0=1,y1=1.1, lwd=1.5, length = 0)

        arrows(x0=final_E_I[length(final_E_I$Inc),1],
               x1=final_E_I[length(final_E_I$Inc),1],
               y0= 1, y1= 1.1,lwd=1.5, length = 0)

        Inc_l <- round(final_E_I[1,1], digits= 1)
        Inc_u <- round(final_E_I[length(final_E_I$Inc),1], digits=1)

        text(x=final_E_I[1,1], y=1, pos=2, Inc_l,cex=1.2)
        text(x=final_E_I[length(final_E_I$Inc),1], y=1, pos=4, Inc_u, cex=1.2)

        #plot text results again if covered by red lines
        if(inc_e_dec[1,1]<40){
          text(x=0, y=3.2,pos=4,text, cex= 0.8)
          text(x=20, y=3.2, pos=4, text2, cex=0.8)
        }
        #recalculate boostrtapped confidence for Edec for plotting confidence margin
        final_E_I_Edec <- final_E_Ibk
        colnames(final_E_I_Edec) <- c("Inc","E","E_dec")
        final_E_I_Edec <- final_E_I_Edec[order(final_E_I_Edec$E_dec),]
        final_E_I_Edec <- final_E_I_Edec[Lconf:Uconf,]

        #plot histogram of E_declination with respect V1 before and after correction
        if(hist==TRUE){
          par(fig=c(0.6,1,0.56,0.99), new=TRUE)
          hist(init_E_Ibk$E_dec, xlim=c(-90,90), breaks= 90,
               axes=FALSE,xlab="",ylab="",col="blue", border="blue", main="")
          par(fig=c(0.6,1,0.56,0.99), new=TRUE)
          hist(final_E_Ibk$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
               breaks = 90, xlab = "", ylab = "",
               main="", cex.axis=0.8,col="red",border ="red")
          abline(v=final_E_I_Edec[1,3],lwd=1, lty=2)
          abline(v=final_E_I_Edec[length(final_E_I_Edec[,3]),3],lwd=1, lty=2)

          #plot lables closer than standard to axes
          title(xlab = "Edec(°)", line=1.9)
          title(ylab = "Frequency", line=1.9)
        }
        #build result file tab
        unstr_boot_stat <- as.data.frame(matrix(ncol = 3,nrow =2))
        rownames(unstr_boot_stat) <- c("Inc", "Edec")
        colnames(unstr_boot_stat) <- c("Mean","Low", "High")
        unstr_boot_stat[1,1] <- inc_nstr
        unstr_boot_stat[1,2] <- final_E_I[1,1]
        unstr_boot_stat[1,3] <- final_E_I[nrow(final_E_I),1]
        unstr_boot_stat[2,1] <- E_nstr
        unstr_boot_stat[2,2] <- final_E_I_Edec[1,3]
        unstr_boot_stat[2,3] <- final_E_I_Edec[nrow(final_E_I_Edec),3]
        return(unstr_boot_stat)
      }

      #performe routine
      Unstrain_Dirs_b <- function(){
        S_vec <- inv_result[[1]]
        unstr_boot_S(unstr_file = unstr_result,nb = input$strboot,S_vec = S_vec,
                     Lin= input$lin,Fol = input$fol,ns = input$incremboot,
                     cross = ifelse(input$brk==2, TRUE, FALSE),
                     EdMAX = ifelse(input$brk==3, TRUE, FALSE),
                     EdMIN = ifelse(input$brk==4, TRUE, FALSE),
                     hist = ifelse(input$strainhist==1, TRUE, FALSE))
      }
      unstr_boot_result <- Unstrain_Dirs_b()
      assign("unstr_boot_result",unstr_boot_result,envir = .GlobalEnv)
      output$unstr_boot_stat <- renderTable({
        unstr_boot_result
      }, rownames = T, digits = 2)
    })

    #perform and plot initial unstrain
    output$unstrainplot <- renderPlot({
      unstrain()
      unstrG <- recordPlot()
      output$unstrbootG <- downloadHandler(
        filename = function() {
          paste(input$fileN_str,"_unstrStat_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 12,height = 11)
          replayPlot(unstrG)
          dev.off()
        }
      )

    }, width = 800,height = 700)
    #perform and plot Unstrain Bootstrap
    output$unstrainboot <- renderPlot({
      unstrain_b()
      unstrBG <- recordPlot()
      output$unstrbootG <- downloadHandler(
        filename = function() {
          paste(input$fileN_strBoot,"_unstrBoot_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 12,height = 11)
          replayPlot(unstrBG)
          dev.off()
        }
      )
      output$unstr_boot_result <- downloadHandler(
        filename = function(){
          paste(input$fileN_strBoot,"unstrStat_",Sys.Date(),".csv",sep = "")
        },
        content = function(file){
          write.csv(unstr_boot_result,file,row.names = T)
        }
      )
    }, width=800, height= 700)
    ############ END OF UNSTRAIN MODULE

    ############ MAGNETIC POLARITY MODULE
    #set reactive file
    Tab_normals <- reactiveValues(list=NULL)

    #main function
    mgstr <- function() {
      #data are always tilt corrected
      DI <- fix_DI(input_file(),coord = 2)

      #fix base and top
      if(is.na(input$baseMS)==FALSE){
        DI <- DI[(DI[,3]>=input$baseMS),]
      }
      if(is.na(input$topMS)==FALSE){
        DI <- DI[(DI[,3]<=input$topMS),]
      }
      #sort position upw
      DI <- DI[order(DI[,3]),]

      #plot magstrat
      if(input$colmgstr==1) colmgstr <- "black"
      if(input$colmgstr==2) colmgstr <- "blue"
      if(input$colmgstr==3) colmgstr <- "green"
      if(input$colmgstr==4) colmgstr <- "pink"
      if(input$colmgstr==5) colmgstr <- "purple"
      if(input$colmgstr==6) colmgstr <- "brown"
      if(input$colmgstr==7) colmgstr <- "red"
      if(input$colmgstr==8) colmgstr <- "yellow"
      if(input$colmgstr==9) colmgstr <- "cyan"
      if(input$colmgstr==10) colmgstr <- "gray"
      if(input$colmgstr==11) colmgstr <- "white"

      #if dec offset is empty then equal 0
      ifelse(is.na(input$Doffset)==T,decoffset <- 0,decoffset <- input$Doffset)

      Tab_normals$list <- PmagDiR::magstrat_DI(DI,lat = input$lat,long = input$long,offset=decoffset, plot_ext = F,POLE = F,E.A. = F,
                                               cex.main = 1.8,cex.lab = 1.6,cex.axis = 1.3,lwd.grid = 1.2,
                                               h_grid=input$hGrid,col = colmgstr,Shiny = T)
    }
    output$magstrat <- renderPlot({
      mgstr()
      #record plot
      mgstrPlot <- recordPlot()
      #Export graphic

      output$mgstr <- downloadHandler(
        filename = function() {
          paste(input$fileN_mgstr,"_mgstr_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, onefile = TRUE,width = 12,height = 10)
          replayPlot(mgstrPlot)
          dev.off()
        }
      )
      #export tab with normal low and high level
      output$revTab <- downloadHandler(
        filename=function() {
          paste(input$fileN_mgstr,"_normal_zone_", Sys.Date(),".csv",sep="")
        },
        content=function(file){
          write.csv(Tab_normals$list,file, row.names = F)
        }
      )
    },width = 900,height = 800)
    ############ END OF MAGNETIC POLARITY MODULE

  }
  shinyApp(ui, server)
}


#service functs
#directions in Cartesian coordinates, and vice versa
s2c <- function(DI,J=1){
  #cart do rad and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI
  data$x <- J*cos(d2r(data[,1]))*cos(d2r(data[,2]))
  data$y <- J*sin(d2r((data[,1])))*cos(d2r((data[,2])))
  data$z <- J*sin(d2r((data[,2])))
  result <- data[,-c(1,2)]
  return(result)
}
c2s <- function(xyz){
  #cart do rad and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- xyz
  data$dec <- r2d(atan2(data[,2],data[,1]))
  data$inc <- r2d(asin(data[,3]/(sqrt((data[,1]^2)+(data[,2]^2)+(data[,3]^2)))))
  result <- data[,-c(1,2,3)]
  return(result)
}

#return the eigenvectors of the AMS inverse matrix for later unstrain
AMS_inv <- function(mat,type="v",prnt=TRUE, Shiny=FALSE){
  library(matlib)
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  AMS <- as.matrix(mat)
  #if type = v AMS is in vector form
  if (type=="v"){
    #AMS in the form of Vi, Vi_dec, Vi_inc
    A <- as.data.frame(matrix(c(AMS[2],AMS[3],
                                AMS[5],AMS[6],
                                AMS[8],AMS[9]),
                              ncol=2,byrow=T))
    colnames(A) <- c("dec","inc")
    A$x <- cos(d2r(A$dec))*cos(d2r(A$inc))
    A$y <- sin(d2r(A$dec))*cos(d2r(A$inc))
    A$z <- sin(d2r(A$inc))
    A_vec <- t(A[,3:5])
    A_val <- c(AMS[1],0,0,0,AMS[4],0,0,0,AMS[7])
    A_val <- matrix(A_val,nrow=3,byrow=T)
    AMS_M <- A_vec%*%A_val%*%inv(A_vec)
  }
  #if type = m AMS is in matrix form
  if (type=="m") {
    #AMS in the form of K11,K22,K33,K12,K23,K13
    AMS_M <- c(AMS[1],AMS[4],AMS[6],
               AMS[4],AMS[2],AMS[5],
               AMS[6],AMS[5],AMS[3])
    AMS_M <- matrix(AMS_M,nrow=3,byrow=T)
  }
  #take full 3x3 matrix
  if(type=="m3x3"){AMS_M <- mat}
  #inverts matrix
  invAMS <- inv(AMS_M)
  inv_AMSe <- eigen(invAMS,symmetric = TRUE)
  AMS_inv <- inv_AMSe$vectors
  AMS_inv_val <- inv_AMSe$values
  AMS_Me <- eigen(AMS_M,symmetric = TRUE)
  AMS_Mval <- AMS_Me$values
  #original anisotropy parameter
  L <- AMS_Mval[1]/AMS_Mval[2]
  F <- AMS_Mval[2]/AMS_Mval[3]
  P <- AMS_Mval[1]/AMS_Mval[3]
  #inverted anisotropy parameter
  Li <- AMS_inv_val[1]/AMS_inv_val[2]
  Fi <- AMS_inv_val[2]/AMS_inv_val[3]
  Pi <- AMS_inv_val[1]/AMS_inv_val[3]
  #calculate inverted direction of axese if shiny==true
  if(Shiny==TRUE){
    V1_inc <- r2d(asin(AMS_inv[3,1]/(sqrt((AMS_inv[1,1]^2)+(AMS_inv[2,1]^2)+(AMS_inv[3,1]^2)))))
    V1_dec <- (r2d(atan2(AMS_inv[2,1],AMS_inv[1,1])))%%360

    V2_inc <- r2d(asin(AMS_inv[3,2]/(sqrt((AMS_inv[1,2]^2)+(AMS_inv[2,2]^2)+(AMS_inv[3,2]^2)))))
    V2_dec <- (r2d(atan2(AMS_inv[2,2],AMS_inv[1,2])))%%360

    V3_inc <- r2d(asin(AMS_inv[3,3]/(sqrt((AMS_inv[1,3]^2)+(AMS_inv[2,3]^2)+(AMS_inv[3,3]^2)))))
    V3_dec <- (r2d(atan2(AMS_inv[2,3],AMS_inv[1,3])))%%360
    AMS_inv_eigen_tab <- round(matrix(c(V1_dec,V2_dec,V3_dec,Li,Fi,V1_inc,V2_inc,V3_inc),ncol = 5,byrow = T),digits=2)
    colnames(AMS_inv_eigen_tab) <- c("V1","V2","V3","L_inv","F_inv")
    rownames(AMS_inv_eigen_tab) <- c("Dec", "Inc")
    AMS_inv_eigen_tab[2,4:5] <- c("","")
  }

  if(prnt==TRUE){
    #print anisotropy parameters
    cat(paste("Anisotropy parameter:
L:",round(L,digits = 4),"
F:", round(F,digits = 4),"
P:", round(P,digits = 4),"
"))
  }

  #returns inverted anisotropy directions if Shiny is FALSE
  if(Shiny==FALSE){return(AMS_inv)}
  #return result list if Shiny is TRUE
  if(Shiny==TRUE){
    result <- list()
    result[[1]] <- AMS_inv
    result[[2]] <- AMS_inv_eigen_tab
    return(result)
  }
}

#Function that rotate geographic dec_inc pair(DI) into bedding coordinates
#if bedding ins not in the file as 3 and 4 column, can be added in function
bed_DI <- function(DI,in_file=TRUE, bed_az,bed_plunge,export=FALSE){
  #functions degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions spherical (Dec=x, Inc=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  data <- DI
  data <- na.omit(data)
  #colnames change if BdecInc are in file or not
  if(in_file==FALSE){
    if(ncol(data)==3) {
      colnames(data) <- c("dec","inc","position")
      pos_exist=TRUE
      }
    if(ncol(data)==2) {
      colnames(data) <- c("dec","inc")
      pos_exist=FALSE
    }
  } else {
    if(ncol(data)==5) {
      colnames(data) <- c("dec","inc","baz","bplunge","position")
      pos_exist=TRUE
      }
    if(ncol(data)==4) {
      colnames(data) <- c("dec","inc","baz","bplunge")
      pos_exist=FALSE
      }
    }

  #sines and cosines if BdecInc are not in file
  if(in_file==FALSE){
    sbd <- -sin(d2r(bed_az))
    cbd <- cos(d2r(bed_az))
    sbi <- sin(d2r(bed_plunge))
    cbi <- cos(d2r(bed_plunge))
  }

  newDI <- as.data.frame(matrix(ncol=2,nrow=0))
  for(i in 1:nrow(data)){
    newDI_p <- as.data.frame(matrix(ncol=2,nrow=1))
    if(in_file==TRUE){
      sbd <- -sin(d2r(data[i,3]))
      cbd <- cos(d2r(data[i,3]))
      sbi <- sin(d2r(data[i,4]))
      cbi <- cos(d2r(data[i,4]))
    }
    x <- s2cx(data[i,1],data[i,2])
    y <- s2cy(data[i,1],data[i,2])
    z <- s2cz(data[i,2])
    xn <- x*(sbd^2+cbd^2*cbi)+
      y*(cbd*sbd*(1-cbi))+
      z*sbi*cbd
    yn <- x*cbd*sbd*(1-cbi)+
      y*(cbd^2+sbd*sbd*cbi)-
      z*sbd*sbi
    zn <- -(x*cbd*sbi-
              y*sbi*sbd-
              z*cbi)
    newdec <- r2d(atan2(yn,xn))
    newdec <- ifelse(newdec<0,newdec+360,newdec)
    newinc <- r2d(asin(zn))
    newDI_p[1,1:2] <- c(newdec,newinc)
    newDI <- rbind(newDI,newDI_p)
  }
  if(pos_exist==TRUE){
    newDI <- cbind(newDI,data$position)
    colnames(newDI) <- c("TCdec","TCinc","position")
  }else{colnames(newDI) <- c("TCdec","TCinc")}
  if(export==TRUE){write.csv(newDI,"tilt_corrected_directions.csv",row.names = FALSE)}
  return(newDI)
}

#check_bipolarity
bip_check <- function(DI){
  #fucnctions deg to rads and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")

  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))

  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))

  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)

  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$values

  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- V1dec%%360

  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  m1ind <- as.numeric(which(data$diff<=90), arr.ind = TRUE)
  m2ind <- as.numeric(which(data$diff>90), arr.ind = TRUE)
  t <- ifelse(length(m2ind)>0 && length(m1ind)>0, TRUE, FALSE)
  return(t)
}

#function that generates resampled Data dec_inc
boots_DI <- function(DI) {
  library("tidyverse", warn.conflicts = FALSE)
  n <- nrow(DI)
  newDI <- DI[sample(n,n,replace = T),]
  return(newDI)
}

#interpolate great circle through directions and return dec inc of pole
circle_DI <- function(DI){
  #degrees to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #coordinate of V3
  V3inc <- r2d(asin(T_vec[3,3]/(sqrt((T_vec[1,3]^2)+(T_vec[2,3]^2)+(T_vec[3,3]^2)))))
  V3dec <- (r2d(atan2(T_vec[2,3],T_vec[1,3])))%%360
  if(V3inc<0){
    V3dec <- V3dec+180
    V3inc <- abs(V3inc)
  }
  #calculate MAD following Kirschvink 1980
  MAD_C <- r2d(atan(sqrt((T_val[3]/T_val[2])+(T_val[3]/T_val[1]))))

  return(c(V3dec,V3inc,MAD_C))
}

#flips all data toward common polarity
common_DI <- function(DI,down=TRUE, export=FALSE,name="common_dirs") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)
  #force V1 down if down=TRUE or up if FALSE
  if(down==TRUE){
    if(V1inc<0){
      V1inc <- -V1inc
      V1dec <- ifelse((V1dec+180)>360,V1dec-180,V1dec+180)
    }
  }
  if(down==FALSE){
    if(V1inc>=0){
      V1inc <- -V1inc
      V1dec <- ifelse((V1dec+180)>360,V1dec-180,V1dec+180)
    }
  }
  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #flips all directions
  data$dec_N <- ifelse(data$diff>90,ifelse((data$dec+180)>360,data$dec-180,data$dec+180),data$dec)
  data$inc_N <- ifelse(data$diff>90,-data$inc,data$inc)
  new_dec_inc <- subset(data,select=c(dec_N,inc_N))
  colnames(new_dec_inc) <- c("dec", "inc")
  new_dec_inc
  #if export==TRUE export flipped data into file
  if(export==TRUE){
    write.csv(new_dec_inc,paste(name, ".csv"), row.names = FALSE)
    # cat(paste("File saved as",name,".csv"))
  }
  return(new_dec_inc)
}

#combine directions and great circle, takes two files with dec inc of directions and poles, return directions on GC
comb_GC_dirs <- function(dirs,poles){
  #cart do rad and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #directions in Cartesian coordinates, return 3 columns and n rowsuni
  s2c <- function(DI){
    #cart do rad and vice versa
    d2r <- function(x) {x*(pi/180)}
    r2d <- function(x) {x*(180/pi)}
    x <- cos(d2r(DI[,1]))*cos(d2r(DI[,2]))
    y <- sin(d2r((DI[,1])))*cos(d2r((DI[,2])))
    z <- sin(d2r((DI[,2])))
    result <- data.frame(matrix(t(c(x,y,z)),ncol = 3))
    colnames(result) <- c("x","y","z")
    return(result)
  }
  c2s <- function(xyz){
    #cart do rad and vice versa
    d2r <- function(x) {x*(pi/180)}
    r2d <- function(x) {x*(180/pi)}
    dec <- r2d(atan2(xyz[1,2],xyz[1,1]))
    inc <- r2d(asin(xyz[1,3]/(sqrt((xyz[1,1]^2)+(xyz[1,2]^2)+(xyz[1,3]^2)))))
    result <- data.frame(matrix(t(c(dec,inc)),ncol = 2))
    colnames(result) <- c("dec","inc")
    return(result)
  }

  #save number of directions
  dirs_num <- nrow(dirs)

  #set initial guess if there are no dirs
  if(dirs_num==0){
    dirs <- data.frame(t(c(sample(0:359,1),sample(-90:90,1))))
  }
  colnames(dirs) <- c("dec","inc")

  #create directions on great circle empty file
  dirs_GC <- data.frame(matrix(ncol = 2,nrow = 0))
  colnames(dirs_GC) <- c("dec","inc")
  n <- 1
  #first_loop
  repeat{
    dirs_f <- PmagDiR::fisher(dirs)
    guess0 <- dirs_f[1,1:2]
    guess0_C <- s2c(guess0)
    #take each pole (row number n) and convert in cart
    pole_C <- s2c(poles[n,1:2])
    #calculate tau
    tau <- (guess0_C[1,1]*pole_C[1,1])+(guess0_C[1,2]*pole_C[1,2])+(guess0_C[1,3]*pole_C[1,3])
    #calculate ro
    ro <- (sqrt(1-(tau^2)))
    #calculate coordinates of dir on great circle
    Xg <- (guess0_C[1,1]-(tau*pole_C[1,1]))/ro
    Yg <- (guess0_C[1,2]-(tau*pole_C[1,2]))/ro
    Zg <- (guess0_C[1,3]-(tau*pole_C[1,3]))/ro
    XYZg <- (matrix(t(c(Xg,Yg,Zg)),ncol=3))
    dir_on_GC <- c2s(XYZg)
    dirs <- rbind(dirs,dir_on_GC)
    dirs_GC <- rbind(dirs_GC,dir_on_GC)
    if(n==nrow(poles)) break
    n <- n+1
  }
  #eliminate first random guess if no dirs were present
  if(dirs_num==0) {dirs <- dirs[-1,]}
  #create new loop where directions from GC are recalculated one by one
  #double loop compares result file with precious result file, break if are equal
  check <- 1
  repeat{
    #l is the first great circle dirs in the complete file
    l <- dirs_num+1
    m <- 1
    #copy dirs file for comparison
    dirs_pre <- dirs
    repeat{
      dirs_t <- dirs[-l,]
      dirs_t_f <- PmagDiR::fisher(dirs_t)
      guess0 <- dirs_t_f[1,1:2]
      guess0_C <- s2c(guess0)
      #take each pole (row number m) and convert in cart
      pole_C <- s2c(poles[m,1:2])
      #calculate tau
      tau <- (guess0_C[1,1]*pole_C[1,1])+(guess0_C[1,2]*pole_C[1,2])+(guess0_C[1,3]*pole_C[1,3])
      #calculate ro
      ro <- (sqrt(1-(tau^2)))
      #calculate coordinates of dir on great circle
      Xg <- (guess0_C[1,1]-(tau*pole_C[1,1]))/ro
      Yg <- (guess0_C[1,2]-(tau*pole_C[1,2]))/ro
      Zg <- (guess0_C[1,3]-(tau*pole_C[1,3]))/ro
      XYZg <- (matrix(t(c(Xg,Yg,Zg)),ncol=3))
      dir_on_GC <- c2s(XYZg)
      dirs[l,] <- dir_on_GC
      dirs_GC[m,] <- dir_on_GC
      if(m==nrow(poles)) break
      l <- l+1
      m <- m+1
    }
    if(all(round(dirs,digits = 3)==round(dirs_pre,digits = 3))) break
    check <- check+1
  }
  #calculate best R
  all_dirs_f <- PmagDiR::fisher(dirs)
  R <- all_dirs_f[1,5]
  #estimate k
  k <- ((2*dirs_num)+nrow(poles)-2)/(2*(dirs_num+nrow(poles)-R))
  #return only directions on GC in this version
  return(dirs_GC)
}

#find crossing point of two dataset of xy coordinates
curve_cross <- function(a, b) {
  colnames(a) <- c("x","y")
  colnames(b) <- c("x","y")
  curve1 <- approxfun(a$x, a$y, rule = 2)
  curve2 <- approxfun(b$x, b$y, rule = 2)

  inters_x <- uniroot(function(x) curve1(x) - curve2(x),
                      c(min(a$x), max(a$x)))$root

  inters_y <- curve2(inters_x)
  cross <- cbind(inters_x,inters_y)
  return(cross)
}

#Dynamic VANDAMME or VGP(45) cutoff (EI before the filter)
#(Physics of the Earth and Planetary Interiors 85;1994)
cut_DI <- function(DI,VD=TRUE,lat=0,long=0,cutoff=40, geo=FALSE,inc_f=TRUE, export=FALSE, name="cut_dirs",Shiny=F){
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #if coordinates are geographic, takes the four columns to calculate TC directions
  if(geo==TRUE){
    DIAP <- DI
    data <- bed_DI((DIAP))
    data <- data[,1:2]
  }else{data <- DI[,1:2]}

  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #enter longitude an latitude and convert radians
  Slat_d <- lat
  Slong_d <- long
  Slat_r <- d2r(Slat_d)
  Slong_r <- d2r(Slong_d)
  #reiterate Vandamme filter applying incl-flattening correction first
  n <- 0
  repeat{
    #number of reiteration
    n <- n+1
    if(inc_f==TRUE){
      #calculate f factor of distribution and f=1 if it is not flattened
      I_E_Edec_f <- ffind(data,f_inc = 0.005)
      f <- ifelse(is.na(I_E_Edec_f[1,4])==TRUE, 1, I_E_Edec_f[nrow(I_E_Edec_f),4])
    }else{f <- 1}

    #add column with inc unflattened, plus different parameters
    unfl_data <- PmagDiR::unflat_DI(data,f)
    data$inc_u <- unfl_data[,2]
    data$dec_r <- d2r(data$dec)
    data$inc_u_r <- d2r(data$inc_u)
    data$Slat_r <- rep(Slat_r)
    data$Slat_d <- rep(Slat_d)
    data$Slong_r <- rep(Slong_r)
    data$Slong_d <- rep(Slong_d)
    #calculate pole colatitude
    data$p_colat_r <- atan(2/tan(data$inc_u_r))
    data$p_colat_d <- r2d(data$p_colat_r)
    #calculate pole latitude
    data$PLat_r <- asin((sin(data$Slat_r)*cos(data$p_colat_r))+
                          (cos(data$Slat_r)*sin(data$p_colat_r)*cos(data$dec_r)))
    data$Plat_d<- r2d(data$PLat)
    #Longitudinal difference between site and pole
    data$LDist_r <- asin((sin(data$p_colat_r)*sin(data$dec_r))/cos(data$PLat_r))
    data$LDist_d <- r2d(data$LDist_r)
    #calculate longitude
    data$PLong_d <- ifelse(cos(data$p_colat_r)<(sin(data$Slat_r)*sin(data$PLat_r)),
                           data$Slong_d+180-data$LDist_d,
                           data$Slong_d+data$LDist_d)
    #long_lat in cartesian coordinates
    data$x <- cos(d2r(data$PLong_d))*cos(data$PLat)
    data$y <- sin(d2r(data$PLong_d))*cos(data$PLat)
    data$z <- sin(data$PLat)
    #average cartesian coordinates
    X_aver <- mean(data$x)
    Y_aver <- mean(data$y)
    Z_aver <- mean(data$z)
    #sum of all values along the axes
    X_sum <- sum(data$x)
    Y_sum <- sum(data$y)
    Z_sum <- sum(data$z)
    #magnitude of average
    B <- sqrt((X_aver^2)+(Y_aver^2)+(Z_aver^2))

    N <- as.numeric(length(data$dec))
    #calculate paleomagnetic pole
    Long_aver <- r2d(atan2(Y_aver,X_aver)) %% 360
    # #corrects for negative declination
    Lat_aver <- r2d(asin(Z_aver/B))
    #PPole long
    data$Pole_long <- rep(Long_aver)
    #PPole lat
    data$Pole_lat <- rep(Lat_aver)
    data$delta <- abs(data$PLong_d-data$Pole_long)
    #calculate angle between Pole and VGP
    data$diff <- r2d(acos((sin(d2r(data$Plat_d))*sin(d2r(data$Pole_lat)))+
                            (cos(d2r(data$Plat_d))*cos(d2r(data$Pole_lat))*cos(d2r(data$delta)))))
    if(VD==TRUE){
      #vandamme filtering calculation
      ASD <- sqrt(sum(((data$diff)^2)/(N-1)))
      A <- (1.8*ASD)+5
    }else{
      A <- cutoff
    }
    #Determine index of data to cut
    VGPcut_t <- as.numeric(which(data$diff>A), arr.ind = T)
    #if there are directions to cut
    if(length(VGPcut_t)!=0){
      #vandamme cut one direction at the time (more loops)
      if(VD==TRUE){
        max_2_cut <- max(data[VGPcut_t,ncol(data)])
        VGPcut <- as.numeric(which(data$diff==max_2_cut,arr.ind = T))
        #otherwise it cuts all that are >A
      }else{VGPcut <- VGPcut_t}
      #cut working data file
      data <- data[-VGPcut,]
      #cut also lines from DI or DIAP file if coordinates are geographic
      if(geo==TRUE){
        DIAP <- DIAP[-VGPcut,]
        }else{DI <- DI[-VGPcut,]}
      #eliminate calculation columns
      data  <- data[,1:2]
      #if VGPcut_temporary is empty, it ends the loop
    }else{break}
  }
  #export cut file
  if(export==TRUE){
    if(geo==TRUE){write.csv(round(DIAP,digits = 2),paste(name,".csv"),row.names = FALSE)}else
    {write.csv(round(DI,digits = 2),paste(name,".csv"),row.names = FALSE)}
  }
  if(Shiny==F){
    cat(paste("Number of reiteration: ", n,"
"))
  }
  #return file
  ifelse(geo==TRUE,
         return(DIAP),
         return(DI))
}

#function that apply cutoff to VGP (not starting from directions, used in web version)
cut_VGP <- function(VGP,VD=TRUE,cutoff=40){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #create file to manipulate
  dat <- VGP
  colnames(dat) <- c("Long","Lat")
  N <- nrow(dat)
  #calculate average
  Pole <- PmagDiR::fisher(dat)
  dat$Plong <- rep(Pole[1,1])
  dat$Plat <- rep(Pole[1,2])
  #calculate distance pole - VGP
  dat$delta <- abs(dat$Long-dat$Plong)
  dat$diff <- r2d(acos((sin(d2r(dat$Lat))*sin(d2r(dat$Plat)))+
                         (cos(d2r(dat$Lat))*cos(d2r(dat$Plat))*cos(d2r(dat$delta)))))
  #apply cutoff
  if(VD==TRUE){
    #vandamme filtering calculation
    repeat{
      ASD <- sqrt(sum(((dat$diff)^2)/(N-1)))
      A <- (1.8*ASD)+5
      VGPcut_t <- as.numeric(which(dat$diff>A), arr.ind = TRUE)
      if(length(VGPcut_t>0)){
        #find maximum deviated VGP and its index
        max_2_cut <- max(dat[VGPcut_t,ncol(dat)])
        VGPcut <- as.numeric(which(dat$diff==max_2_cut,arr.ind = T))
        dat <- dat[-VGPcut,]
      }else(break)
    }
  }else if(VD==FALSE){
    A <- cutoff
    VGPcut <- as.numeric(which(dat$diff>A), arr.ind = TRUE)
    if(length(VGPcut>0)){dat <- dat[-VGPcut,]}
  }
  VGP <- dat[,1:2]
  return(VGP)
}


#convert VGPs and site latitude and longitude in directions
DI_from_VGP <- function(VGPs, lat, long, export=FALSE,name="Directions") {
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  VGPs <- na.omit(VGPs)
  directions <- as.data.frame(matrix(nrow = 0,ncol = 2))
  for (i in 1:nrow(VGPs)){
    plong <- VGPs[i,1] %% 360
    plat <- VGPs[i,2]
    long <- long %% 360
    signdec <- 1
    Dphi <- abs(plong - long)
    if (Dphi != 0) signdec <- (plong - long) / Dphi
    if (lat == 90) lat <- 89.99
    thetaS <- d2r(90 - lat)
    thetaP <- d2r(90 - plat)
    Dphi <- d2r(Dphi)
    cosp <- cos(thetaS) * cos(thetaP) + sin(thetaS) * sin(thetaP) * cos(Dphi)
    thetaM <- acos(cosp)
    cosd <- (cos(thetaP) - cos(thetaM) * cos(thetaS)) / (sin(thetaM) * sin(thetaS))
    C <- abs(1 - cosd^2)
    dec <- ifelse(C != 0,-atan(cosd / sqrt(abs(C))) + (pi / 2),dec <- acos(cosd))
    if (-pi < signdec * Dphi && signdec < 0) dec <- 2 * pi - dec
    if (signdec * Dphi > pi) dec <- 2 * pi - dec
    dec <- r2d(dec) %% 360
    inc <- r2d(atan2(2 * cos(thetaM), sin(thetaM)))
    DecInc <- as.data.frame(t(c(dec,inc)))
    directions <- rbind(directions,DecInc)
  }
  colnames(directions) <- c("Dec","Inc")
  #export csv with directions if requested
  if(export==TRUE) {
    write.csv(round(directions,digits = 2),paste(name,".csv"),row.names = FALSE)
  }
  return(directions)
}

#function that calculate E_I couples of data and plot bootstrapped statistics
EI_boot <- function(DI,nb=1000,conf=95,export=TRUE, name="EI_boot_plot") {
  data <- DI
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec", "inc")
  Inc_E_real <- inc_E_finder(data)
  Inc_E_real$V1inc <- abs(Inc_E_real$V1inc)
  Inc_E <- as.data.frame(matrix(ncol=3,nrow=0))
  cat(paste("Simulating",nb,"detasets.
"))
  for (i in 1:nb) {
    dataprov <- boots_DI(data)
    I_E_Ed <- inc_E_finder(dataprov)
    I_E_Ed$V1inc <- abs(I_E_Ed$V1inc)
    Inc_E <- rbind(Inc_E, I_E_Ed)
    if(((i%%50)==0)==TRUE) {
      cat(paste(i,"simulations out of",nb,"done
"))
    }
  }
  colnames(Inc_E) <- c("Inc","E","E_dec")
  Inc_E <- Inc_E[order(Inc_E$E),]
  Inc_E_bk <- Inc_E
  Inc_E <- Inc_E_bk

  confn <- conf/100
  num <- round((nb*(1-confn))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  Inc_E <- Inc_E[Lconf:Uconf,]

  N <- as.character(length(data[,1]))
  Inc <- format(round(Inc_E_real$V1inc,1),nsmall=1)
  Ecut <- format(round(Inc_E_real$E,2),nsmall=2)
  V2 <- format(round(Inc_E_real$DV1V2,1),nsmall=1)

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(Inc_E_real$E>3.5, ceiling(Inc_E_real$E),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #Plot Bootstrapped data
  points(x=Inc_E$Inc,
         y=Inc_E$E,
         pch=16,
         col=rgb(0, 0, 1, 0.15),
         cex=0.6)

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #Plot real data E-I couple and values
  points(x=Inc_E_real$V1inc,y=Inc_E_real$E,pch=21,
         col="black", bg="white", cex=1.2)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #plot histogram of Dec top right
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)

  hist(Inc_E$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
       breaks = 90, xlab = "", ylab = "",
       main="", cex.axis=0.8, col= "blue", border="blue")
  #plot labels closer than standard to axes
  title(xlab = "Edec(°)", line=1.9, cex=0.3)
  title(ylab = "Frequency", line=2, cex=0.1)

  #plot confidence margin of declination if between -45 and 45
  if(Inc_E_real$DV1V2>-45 && Inc_E_real$DV1V2<45){
    Inc_E_Edec <- Inc_E_bk
    Inc_E_Edec <- Inc_E_Edec[order(Inc_E_Edec[,3]),]
    Inc_E_Edec <- Inc_E_Edec[Lconf:Uconf,]
    low_dec <- Inc_E_Edec[1,3]
    up_dec <- Inc_E_Edec[length(Inc_E_Edec[,3]),3]
    abline(v=low_dec,lwd=1,lty=2)
    abline(v=up_dec,lwd=1,lty=2)
  }
  #calculate low and high inclination error
  Inc_E <- Inc_E_bk
  Inc_E <- Inc_E[order(Inc_E[,1]),]
  Inc_E <- Inc_E[Lconf:Uconf,]
  low_inc <- Inc_E[1,1]
  up_inc <- Inc_E[length(Inc_E[,1]),1]

  if(bip_check(data)==TRUE){
    #plot equal area data two modes
    par(fig=c(0.55,1,0,0.6), new=TRUE)
    plot_DI(data,title="Original directions")
  }

  #plot equal area single mode
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(data,single_mode = TRUE,title = "Single mode directions")

  if(Inc_E_real$DV1V2>-45 && Inc_E_real$DV1V2<45){
    results <- as.data.frame(matrix(ncol= 9, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec","Low_E_dec","High_E_dec")
    results$Inc <- Inc
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- Ecut
    results$Low_E <- round(min(Inc_E$E),digits= 2)
    results$High_E <- round(max(Inc_E$E),digits= 2)
    results$E_dec <- V2
    results$Low_E_dec <- round(low_dec,digits= 2)
    results$High_E_dec <- round(up_dec,digits=2)
  }else{
    results <- as.data.frame(matrix(ncol= 7, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec")
    results$Inc <- Inc
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- Ecut
    results$Low_E <- round(min(Inc_E$E),digits=2)
    results$High_E <- round(max(Inc_E$E),digits= 2)
    results$E_dec <- V2
  }
  print(results, row.names = FALSE)
  #reset screen
  par(fig=c(0,1,0,1))
  #export if requested
  if(export==TRUE){
    cat("
Results saved as ", paste(name,".csv"),"
Graph saved as", paste(name,".pdf"),"
")
    write.csv(results,file=paste(name,".csv"),row.names = FALSE)
    save_pdf(name=paste(name,".pdf"))
  }
}

#function that calculated DeltaDec and DeltaInc
ellips_DI <- function(DI,lat,long,export=FALSE){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #cut NA rows
  DI <- na.omit(DI)
  #calculate average dir and paleolat in radians
  aver_DI <- fisher(DI)
  aver_inc <- aver_DI[1,2]
  lat_r <- atan((tan(d2r(aver_inc)))/2)
  #calculate vgps
  poles <- VGP_DI(DI,in_file=FALSE,lat=lat,long=long,export=F,type="VGPs",Prnt=FALSE)
  #calculate A95
  PPole <- fisher(poles)
  A95 <- PPole[1,3]
  #calculated ∂Dec and ∂Inc
  D_dec <- r2d(asin((sin(d2r(A95)))/cos(lat_r)))
  D_inc <- r2d((2*d2r(A95))/(1+(3*(sin(lat_r))^2)))
  result <- as.data.frame(matrix(nrow = 1,ncol=4))
  colnames(result) <- c("dec","inc","delta_dec","delta_inc")
  result[1,1] <- aver_DI[1,1]
  result[1,2] <- aver_DI[1,2]
  result[1,3] <- D_dec
  result[1,4] <- D_inc
  result$N <- nrow(DI)
  #esport results if requested
  if(export==TRUE){
    write.csv(round(result,digits = 2),"confidence_ellipse.csv",row.names = FALSE)
    cat("Confidence ellipse:
")
    print(result)
  }
  return(result)
}

#plot bimodal elliptical confidence (calculated from A95 inversion) from dec_inc  and print results on console
ellips_plot <- function(DI,lat=0,long=0, plot=TRUE, on_plot=TRUE, col_d="red",col_u="white",col_l="black",symbol="c", text=FALSE,export=TRUE,save=FALSE,name="ellipse"){
  #degrees to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  data <- DI
  #cut lines with empty cells
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T,symmetric = TRUE)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- V1dec%%360
  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  if(any(data$diff<=90)){
    mode1 <- as.data.frame(data$dec[data$diff<=90])
    mode1$inc <- data$inc[data$diff<=90]
    colnames(mode1) <- c("dec","inc")
  }
  if(any(data$diff>90)){
    mode2 <- as.data.frame(data$dec[data$diff>90])
    mode2$inc <- data$inc[data$diff>90]
    colnames(mode2) <- c("dec","inc")
  }
  #calculate ellipses
  if(exists("mode1")==TRUE) {ellips_M1 <- ellips_DI(mode1, lat=lat, long=long)}
  if(exists("mode2")==TRUE) {ellips_M2 <- ellips_DI(mode2, lat=lat, long=long)}

  if(plot==TRUE){
    if(on_plot==FALSE){
      par(fig=c(0,1,0,1), new= FALSE)
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
           xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
      equalarea()
    }
    if(exists("mode1")==TRUE){generate_ellips(ellips_M1[1,1],ellips_M1[1,2],ellips_M1[1,3],ellips_M1[1,4],
                                              on_plot = TRUE,symbol=symbol, col_d = col_d,
                                              col_u=col_u,col_l=col_l)}
    if(exists("mode2")==TRUE){generate_ellips(ellips_M2[1,1],ellips_M2[1,2],ellips_M2[1,3],ellips_M2[1,4],
                                              on_plot = TRUE,symbol=symbol, col_d = col_d,
                                              col_u=col_u,col_l=col_l)}
  }
  data_M12 <- common_DI(data)
  ellips_M12 <- ellips_DI(data_M12,lat=lat, long=long)
  #plot text with results
  N <- ellips_M12[1,5]
  Dec <- round(ellips_M12[1,1],digits=2)
  Inc <- round(ellips_M12[1,2],digits=2)
  Delta_dec <- round(ellips_M12[1,3],digits=2)
  Delta_inc <- round(ellips_M12[1,4],digits=2)

  if(any(data$diff<=90)) {
    cat("Ellipse Mode 1:
")
    print(round(ellips_M1, digits=2), row.names = FALSE)
    if(export==TRUE){write.csv(round(ellips_M1, digits=2),paste(name,"_mode_1.csv"), row.names = FALSE)}
  }
  if(any(data$diff>90)) {
    cat("Ellipse Mode 2:
")
    print(round(ellips_M2,digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(ellips_M2,digits=2)),paste(name,"_mode_2.csv"), row.names = FALSE)}
  }
  if(any(data$diff>90)) {
    cat("Ellipse common mode:
")
    print(round(ellips_M12, digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(ellips_M12, digits=2)),paste(name,"_mode_1&2.csv"), row.names = FALSE)}
  }
  if (text==TRUE){
    #plot text if true
    par(fig=c(0,1,0,1), new=TRUE)
    text <- paste("N: ",N,"
Dec: ", Dec,"
Inc: ", Inc,"
∆ dec: ", Delta_dec,"
∆ inc: ", Delta_inc)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

    text(x=0.79, y=0.02,pos=4,text, cex= 0.85)
    cat("\nDo not attempt to plot other directions or Fisher mean on the same diagram if text option is set TRUE.\n")
  }

  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function plotting equal area net
equalarea <- function(title="") {
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #plot external circle
  frame_dec = 0:360
  frame_inc=rep(0,length(frame_dec))
  x = a2cx(frame_inc,frame_dec)
  y = a2cy(frame_inc,frame_dec)
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  lines(x, y, col = "black")
  #plot 30deg and 60deg circles
  tdegD <- 0:360
  tdegI <- rep(30, length(tdegD))
  sdegI <- rep(60, length(tdegD))
  tdegx <- a2cx(tdegI,tdegD)
  tdegy <- a2cy(tdegI,tdegD)
  sdegx <- a2cx(sdegI,tdegD)
  sdegy <- a2cy(sdegI,tdegD)
  points(x=tdegx,y=tdegy,type="l", col="gray")
  points(x=sdegx,y=sdegy,type="l", col="gray")
  #plot lines every 45deg
  xd1_1 <- a2cx(0,225)
  xd1_2 <- a2cx(0,45)
  yd1_1 <- a2cy(0,225)
  yd1_2 <- a2cy(0,45)
  xd2_1 <- a2cx(0,-45)
  xd2_2 <- a2cx(0,135)
  yd2_1 <- a2cy(0,-45)
  yd2_2 <- a2cy(0,135)
  lines(c(-1, 1), c(0, 0), col = "gray")
  lines(c(0, 0), c(-1, 1), col = "gray")
  lines(c(xd1_1, xd1_2), c(yd1_1, yd1_2), col = "gray")
  lines(c(xd2_1, xd2_2), c(yd2_1, yd2_2), col = "gray")
  title(xlab = title, line=0.2, cex=0.1)
}

#Function that correct inclination shallowing after tk03.GAD model
ffind_boot <- function(DI,confidence=95,nb=1000, f_increment=0.01,export=TRUE,return=TRUE, name="Unflattened_dirs") {
  data <- DI[,1:2]
  data <- na.omit(data)
  N <- nrow(data)
  colnames(data) <- c("dec", "inc")
  #calculate E-I of real data
  Inc_E_R <- inc_E_finder(data)
  Inc_E_R$V1inc <- abs(Inc_E_R$V1inc)
  Inc <- round(Inc_E_R$V1inc, digits=1)
  Ecut <- round(Inc_E_R$E, digits=2)
  Edec <- round(Inc_E_R$DV1V2, digits=1)
  cat("Calculating precise inclination flattening of raw data.

")

  #calculate E-I correction sequence of real data
  Seq_I_E_R <- ffind(data,f_inc = 0.0005)
  colnames(Seq_I_E_R) <- c("V1inc","E","DV1V2","f")
  alert <- ifelse(length(Seq_I_E_R$V1inc)==1,"y","n")
  Ffinal <- round(Seq_I_E_R[length(Seq_I_E_R$f),4], digits=2)
  Inc_f <- round(Seq_I_E_R[length(Seq_I_E_R$V1inc),1], digits=1)
  Efinal <- round(Seq_I_E_R[length(Seq_I_E_R$E),2], digits=2)
  Edec_f <- round(Seq_I_E_R[length(Seq_I_E_R$DV1V2),3], digits=1)
  f <- min(Seq_I_E_R$f)
  if(alert=="y") f <- 1
  unf_data <- unflat_DI(data,f)
  colnames(unf_data) <- c("dec","inc")

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  plot(NA, xlim= c(0,90), ylim= c(1,3.5), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #Plot real data E-I and correction of real data
  points(x=Seq_I_E_R$V1inc, y= Seq_I_E_R$E, type= "l", col="red", lwd=3)
  points(x=Inc_E_R$V1inc,y=Inc_E_R$E,pch=21,
         col="black", bg="blue", cex=1.2)
  points(x=Inc_f,y=Efinal,
         pch=21, col="black", bg="red", cex=1.2)

  text <- paste("N:", N, "
Inc:", Inc,"
E:", Ecut,"
Edec:",Edec)

  text2 <- paste("f:", Ffinal, "
Inc_Unfl:", Inc_f, "
E_Unfl:", Efinal, "
Edec_Unfl:", Edec_f)

  text(x=0, y=3.2,pos=4,text, cex= 0.8)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)
  if (alert=="y"){
    text3 <- "Distribution not flattened"
    text(x=0, y=3, pos=4,text3,cex=1)
  }

  #plot equal area data
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(data,single_mode = TRUE,title = "Original directions")

  #plot equal area unflat data
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(unf_data,single_mode = TRUE, col_d ="red",col_u = "pink",
          title="Unflattened directions")

  #create files for initial and final readings for the histograms
  init_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  final_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  colnames(init_E_I) <- c("Inc","E","E_dec")
  colnames(final_E_I)<- c("Inc","E","E_dec")
  if(alert=="y") par(fig=c(0,1,0,1))
  if(alert=="y") stop("
DISTRIBUTION NOT FLATTENED.")


  cat("Bootstrapping.
Simulation ends when", nb, "valid pseudosamples are saved.

")
  n <- 0
  par(fig=c(0,0.7,0,1), new=TRUE)
  plot(NA, xlim= c(0,90), ylim= c(1,3.5), xaxt="n",yaxt="n",
       xlab="", ylab="", axes=FALSE)
  repeat {
    n <- n+1
    Seq_I_E_B <- as.data.frame(matrix(ncol=3,nrow=0))
    dataprov <- boots_DI(data)
    Seq_I_E_B <- ffind(dataprov, f_inc = f_increment)
    #plot bootstrapped lines
    points(x=Seq_I_E_B$V1inc, y= Seq_I_E_B$E,
           type= "l", col=rgb(1, 0, 0, 0.15), lwd=1)
    i_E_I <- Seq_I_E_B[1,]
    f_E_I <- Seq_I_E_B[nrow(Seq_I_E_B),]
    colnames(i_E_I) <- c("Inc","E","E_dec")
    colnames(f_E_I)<- c("Inc","E","E_dec")

    #isolate initial and final readings for histograms
    init_E_I <- rbind(init_E_I,i_E_I)
    final_E_I <- rbind(final_E_I,f_E_I)
    init_E_I <- na.omit(init_E_I)
    final_E_I <- na.omit(final_E_I)
    if(((n%%50)==0)==TRUE) {
      cat(paste(n,"simulations done and",(nrow(final_E_I)),"pseudosamples saved
"))
    }

    if(nrow(final_E_I)==nb) {
      cat(paste("Saved",(length(final_E_I[,1])), "pseudosamples after", n,"simulations
"))
      break
    }
  }
  #replot real data with different color
  points(x=Seq_I_E_R$V1inc, y= Seq_I_E_R$E, type= "l", col="yellow", lwd=3)
  points(x=Inc_E_R$V1inc,y=Inc_E_R$E,pch=21,
         col="black", bg="blue", cex=1.2)
  points(x=x, y= y, type= "l", col="blue", lwd=3)
  points(x=Inc_f,y=Efinal,
         pch=21, col="black", bg="red", cex=1.2)


  #replot results in case covered by boostrapps
  text(x=0, y=3.2,pos=4,text, cex= 0.8)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  colnames(init_E_I) <- c("Inc","E","E_dec","f")
  colnames(final_E_I)<- c("Inc","E","E_dec","f")
  final_E_I <- final_E_I[order(final_E_I$Inc),] #order final results by inclination
  final_E_Ibk <- final_E_I

  conf <- confidence/100
  num <- round((nb*(1-conf))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  final_E_I <- final_E_I[Lconf:Uconf,]    #cut bootstrapped results for 95% confidence

  #draw two lines for 95% confidence margin
  arrows(x0=final_E_I[1,1],x1=final_E_I[1,1],
         y0=1,y1=final_E_I[1,2], length = 0,lty=2)

  arrows(x0=final_E_I[length(final_E_I$Inc),1],
         x1=final_E_I[length(final_E_I$Inc),1],
         y0= 1, y1= final_E_I[length(final_E_I$Inc),2],
         length = 0,lty=2)
  Inc_l95 <- round(final_E_I[1,1], digits= 1)
  Inc_u95 <- round(final_E_I[length(final_E_I$Inc),1], digits=1)

  text(x=final_E_I[1,1], y=1, pos=2, Inc_l95)
  text(x=final_E_I[length(final_E_I$Inc),1], y=1, pos=4, Inc_u95)


  #plot histogram of E_declination with respect V1 before and after correction
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
  hist(init_E_I$E_dec, xlim=c(-90,90), breaks= 90,
       axes=FALSE,xlab="",ylab="",col="blue", border="blue", main="")
  #plot lables closer than standard to axes
  title(xlab = "Edec(°)", line=1.9, cex=0.2)
  title(ylab = "Frequency", line=1.9,cex=0.2)
  #after
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
  hist(final_E_Ibk$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
       breaks = 90, xlab = "", ylab = "",
       main="", cex.axis=0.8,col="red",border ="red")

  #recalculate boostrtapped confidence for Edec for plotting confidence margin
  final_E_I_Edec <- final_E_Ibk
  colnames(final_E_I_Edec) <- c("Inc","E","E_dec")
  final_E_I_Edec <- final_E_I_Edec[order(final_E_I_Edec$E_dec),]
  final_E_I_Edec <- final_E_I_Edec[Lconf:Uconf,]
  abline(v=final_E_I_Edec[1,3],lwd=1, lty=2)
  abline(v=final_E_I_Edec[length(final_E_I_Edec[,3]),3],lwd=1, lty=2)

  #plot histogram of inclination before and after correction
  par(fig=c(0.35,0.69,0.25,0.68), new=TRUE)
  hist(init_E_I$Inc,xlim=c(0,90), breaks=90,
       axes=FALSE,xlab="",ylab="",col="blue",border="blue", main="")
  par(fig=c(0.35,0.69,0.25,0.68), new=TRUE)
  hist(final_E_Ibk$Inc, xlim=c(0,90), xaxp=c(0,90,6),
       breaks=90,xlab = "", ylab = "",
       main="",cex.axis=0.8,col="red",border ="red")
  abline(v=final_E_I[1,1],lwd=1,lty=2)
  abline(v=final_E_I[length(final_E_I[,1]),1],lwd=1,lty=2)
  title(xlab = "Inc(°)", line=1.9, cex=0.2)
  title(ylab = "Frequency", line=1.9, cex=0.2)

  stat <- as.data.frame(matrix(ncol=1,nrow=1))
  colnames(stat) <- c("N")
  stat$N <- N
  stat$Inc <- Inc
  stat$E <- Ecut
  stat$Edec <- Edec
  stat$f <- round(f,digits=2)
  stat$Inc_unfl <- Inc_f
  stat$Low_inc <- Inc_l95
  stat$Hign_inc <- Inc_u95
  stat$E_unfl <- Efinal
  stat$Edec_unfl <- Edec_f
  stat$Edec_low <- round(final_E_I_Edec[1,3],digits = 1)
  stat$Edec_high <- round(final_E_I_Edec[length(final_E_I_Edec[,3]),3],digits = 1)

  par(fig=c(0,1,0,1))
  if(export==TRUE){
    cat("
Unflattened directions and statistics saved as .csv file
Graph saved as",paste(name,".pdf"),"

")
    write.csv(round(unf_data,digits = 2),paste(name,".csv"), row.names = FALSE)
    write.csv(stat,paste(name,"_statistic.csv"), row.names=FALSE)
    save_pdf(name=paste(name,".pdf"))
  }
  if(return==TRUE){return(unf_data)}
}

#flattening factor finder function from Dec Inc, results in Inc, E, and E declination
ffind <-function(DI, f_inc=0.005) {
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  inc_E_seq <- as.data.frame(matrix(ncol=4,nrow=0))
  f <- 1.00
  #loop that stops only when tk03 line is crossed from below
  repeat {
    f <- f-f_inc
    data_unfl <- unflat_DI(data,f)
    inc_E_prov <- inc_E_finder(data_unfl)
    inc_E_prov$V1inc <- abs(inc_E_prov$V1inc)
    inc_E_prov$f <- f
    inc_E_seq <- rbind(inc_E_seq,inc_E_prov)
    if(nrow(inc_E_seq)>1){
      E_lim_low <- round(tk03(inc_E_seq[nrow(inc_E_seq)-1,1]), digits = 6)
      E_lim_high <- round(tk03(inc_E_seq[nrow(inc_E_seq),1]),digits = 6)
      if(round(max(inc_E_seq[nrow(inc_E_seq),1],digit=1)>89.5)) break
      if(inc_E_seq[nrow(inc_E_seq)-1,2]< E_lim_low && inc_E_seq[nrow(inc_E_seq),2]>=E_lim_high) break
    }
  }
  #next return data only when E goes below tk03 line
  Emin <- min(inc_E_seq$E)
  Imin <- inc_E_seq$V1inc[inc_E_seq$E==min(inc_E_seq$E)]
  Eminlim <- tk03(Imin)
  if(round(max(inc_E_seq$V1inc),digit=1)>89.5) {return(as.data.frame(t(c(NA, NA, NA, NA))))}
  if(Emin>Eminlim) {return(as.data.frame(t(c(NA, NA, NA, NA))))} else {return(inc_E_seq)}
}

#plot bimodal fisher from dec_inc and print results on console
fisher_plot <- function(DI, plot=TRUE, on_plot=TRUE,col_d="red",col_u="white",col_l="black",symbol="c",text=FALSE,export=TRUE,save=FALSE,name="Fisher_mean") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T,symmetric = TRUE)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)
  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  if(any(data$diff<=90)){
    mode1 <- as.data.frame(data$dec[data$diff<=90])
    mode1$inc <- data$inc[data$diff<=90]
    colnames(mode1) <- c("dec","inc")
  }
  if(any(data$diff>90)){
    mode2 <- as.data.frame(data$dec[data$diff>90])
    mode2$inc <- data$inc[data$diff>90]
    colnames(mode2) <- c("dec","inc")
  }
  if(exists("mode1")==TRUE) {fisher_M1 <- fisher(mode1)}
  if(exists("mode2")==TRUE) {fisher_M2 <- fisher(mode2)}
  if(plot==TRUE){
    if(on_plot==FALSE){
      par(fig=c(0,1,0,1), new= FALSE)
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
           xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
      equalarea()
      }
    if(exists("mode1")==TRUE){plot_a95(fisher_M1[1,1],fisher_M1[1,2],fisher_M1[1,3],
                                       on_plot = TRUE,symbol=symbol, col_d = col_d,
                                       col_u=col_u,col_l=col_l)}
    if(exists("mode2")==TRUE){plot_a95(fisher_M2[1,1],fisher_M2[1,2],fisher_M2[1,3],
                                       on_plot = TRUE,symbol=symbol, col_d = col_d,
                                       col_u=col_u,col_l=col_l)}
  }
  data_M12 <- common_DI(data)
  fisher_M12 <- fisher(data_M12)
  #plot text with results
  Dec <- round(fisher_M12[1,1],digits=2)
  Inc <- round(fisher_M12[1,2],digits=2)
  a <- round(fisher_M12[1,3],digits=2)
  N <- round(fisher_M12[1,4],digits=2)

  if(any(data$diff<=90)) {
    cat("fisher Mode 1:
")
print(round(fisher_M1, digits=2), row.names = FALSE)
if(export==TRUE){write.csv(round(fisher_M1, digits=2),paste(name,"_mode_1.csv"), row.names = FALSE)}
  }
  if(any(data$diff>90)) {
    cat("fisher Mode 2:
")
    print(round(fisher_M2,digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(fisher_M2,digits=2)),paste(name,"_mode_2.csv"), row.names = FALSE)}
  }
  if(exists("fisher_M1")==TRUE | exists("fisher_M2")==TRUE) {
    cat("fisher common mode:
")
    print(round(fisher_M12, digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(fisher_M12, digits=2)),paste(name,"_mode_1&2.csv"), row.names = FALSE)}
  }

  #plot text in figure if requested
  if (text==TRUE){
    #plot text if true
    par(fig=c(0,1,0,1), new=TRUE)
    text <- paste("N: ",N,"
Dec: ", Dec,"
Inc: ", Inc,"
a95%: ", a)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

    text(x=0.79, y=0.02,pos=4,text, cex= 0.85)
    cat("\nDo not attempt to plot other directions or Fisher mean on the same diagram if text option is set TRUE.\n")
  }

  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function that return fisher statistic from dec_inc
fisher <- function(DI, export=FALSE, name="fisher_mean"){
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #add various columns to data, with conversions to rads, lat and long
  data$dec_r <- d2r(data$dec)
  data$inc_r <- d2r(data$inc)
  #long_lat in cartesian coordinates
  data$x <- cos(data$dec_r)*cos(data$inc_r)
  data$y <- sin(data$dec_r)*cos(data$inc_r)
  data$z <- sin(data$inc_r)
  #average cartesian coordinates
  X_aver <- mean(data$x)
  Y_aver <- mean(data$y)
  Z_aver <- mean(data$z)
  #sum of all values along the axes
  X_sum <- sum(data$x)
  Y_sum <- sum(data$y)
  Z_sum <- sum(data$z)
  #magnitude of average
  B <- sqrt((X_aver^2)+(Y_aver^2)+(Z_aver^2))
  #Fisher (1953) parameters R,K, a95
  N <- length(data$dec)
  R <- sqrt(X_sum^2+Y_sum^2+Z_sum^2)
  K <- (N-1)/(N-R)
  a95 <- r2d(acos(1-(((N-R)/R)*(((1/0.05)^(1/(N-1)))-1))))
  Dec_aver <- r2d(atan2(Y_aver,X_aver))
  #corrects for negative declination
  Dec_aver <- ifelse(Dec_aver<0,Dec_aver+360,Dec_aver)
  Inc_aver <- r2d(asin(Z_aver/B))
  result <- as.data.frame(matrix(ncol=6,nrow=1))
  colnames(result) <- c("dec","inc","a95","N","R","k")
  result[1,1:6] <- c(Dec_aver,Inc_aver,a95,N,R,K)
  if(export==TRUE){write.csv(round(result,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(result)
}

#function generating fisher distributed data
fisher_generator <- function(N,k,lon,lat,k_tol){
  #sub-function generating random long lat
  fisherDiR <- function(k){
    r2d <- function(x) {x*(180/pi)}
    L <- exp(-2*k)
    a <- runif(1)*(1-L)+L
    f <- sqrt(-log(a)/(2*k))
    latitude <- 90-r2d(2*asin(f))
    longitude <- r2d(2*pi*runif(1))
    return(c(longitude, latitude))
  }
  #reiterate until k is within the tolerance
  repeat{
    result <- data.frame(matrix(ncol = 2,nrow = 0))
    colnames(result) <- c("lon", "lat")
    for(i in 1:N){
      decinc_temp <- data.frame(t(fisherDiR(k)))
      colnames(decinc_temp) <- c("lon", "lat")
      result <- rbind(result,decinc_temp)
    }
    AverLongLat <- PmagDiR::fisher(result)
    fixed_data <- PmagDiR::bed_DI(result,in_file = F,
                                  bed_az = (AverLongLat[1,1]+180)%%360,bed_plunge =90-AverLongLat[1,2])
    final_VGP <- PmagDiR::bed_DI(fixed_data,in_file = F,
                                 bed_az = lon,bed_plunge =90-lat)
    colnames(final_VGP) <- c("Long","Lat")
    stat <- PmagDiR::fisher(final_VGP)
    k_test <- stat[1,6]
    #check for tolerance
    if(k_tol==0) break
    if(abs(k_test-k)<=k_tol) break
  }
  return(final_VGP)
}

#flat directions with given f
flat_DI <- function(DI,f=1,export=FALSE,name="flattened_dirs") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  DI[,2]<- r2d(atan(tan(d2r(DI[,2]))*f))
  if(export==TRUE){write.csv(round(DI,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(DI)
}

#flip directions to antipodal
flip_DI <- function(DI,export=FALSE,name="flipped_dirs"){
  dat <- DI[,1:2]
  dat <- na.omit(dat)
  colnames(dat) <- c("dec","inc")
  dat_fl <- dat
  dat_fl$dec <- ifelse((dat$dec+180)>360,dat$dec-180,dat$dec+180)
  dat_fl$inc <- -dat$inc
  if(export==TRUE){write.csv(round(dat_fl,digits=2),paste(name, ".csv"),row.names = FALSE)}
  return(dat_fl)
}

#function that plots points on a KavrayskiyVII geographic map
geo_point <- function(S_file=FALSE,symbol="c",col="red",center=0,grid=30,A95=FALSE,fill_A=TRUE,export=TRUE,on_plot=FALSE){
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy in KavrayskiyVII projection
  c2x <- function(lon,lat) {((3*d2r(lon))/2)*(sqrt((1/3)-((d2r(lat)/pi)^2)))}
  c2y <- function(lat) {d2r(lat)}

  #functions spherical (lon=x, lat=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}

  #function that draw circle around point
  circle <- function(longitude,latitude,A,fill=FALSE){
    #empty new circle file
    circle <- as.data.frame(matrix(ncol=2,nrow=0))
    #loop that create A95 and rotate it around new coordinate (dec, inc)
    for (i in seq(0,360,0.5)){
      #temporary circle point
      circleP <- as.data.frame(matrix(ncol=2,nrow=1))
      S_lon <- longitude
      S_lat <- latitude
      A <- A

      #trigonometric parameters for rotation
      sbd <- -sin(d2r(S_lon))
      cbd <- cos(d2r(S_lon))
      sbi <- sin(d2r(90-S_lat))
      cbi <- cos(d2r(90-S_lat))

      # cartesian coordinates of the confidence circle
      x <- s2cx(i,(90-A))
      y <- s2cy(i,(90-A))
      z <- s2cz(90-A)

      #new rotated coordinate
      xn <- x*(sbd^2+cbd^2*cbi)+
        y*(cbd*sbd*(1-cbi))+
        z*sbi*cbd
      yn <- x*cbd*sbd*(1-cbi)+
        y*(cbd^2+sbd*sbd*cbi)-
        z*sbd*sbi
      zn <- -(x*cbd*sbi-
                y*sbi*sbd-
                z*cbi)
      #converted to spherical
      newlon <- r2d(atan2(yn,xn))
      newlon <- ifelse(newlon>180,newlon-360,newlon)
      newlon <- ifelse(newlon<(-180),newlon+360,newlon)
      newlat <- r2d(asin(zn))
      circleP[1,1:2] <- c(newlon,newlat)
      circle <- rbind(circle,circleP)
    }
    colnames(circle) <- c("lon","lat")

    #next is to divide circles into two polygons if crosses the end of map
    circle1 <- as.data.frame(matrix(nrow = 0,ncol = 2))
    colnames(circle1) <- c("lon","lat")
    #creates the breaking line
    breaker <- as.data.frame(t(c(NA,NA)))
    colnames(breaker) <- c("lon","lat")
    #when two points are on the different side of the map, based on 355° distance, it put break in between
    for(l in 2:nrow(circle)){
      if(abs(circle[l-1,1]-circle[l,1])>355){
        provv1 <- circle[l-1,]
        provv2 <- circle[l,]
        circle1 <- rbind(circle1,provv1,breaker,provv2)
      }else{circle1 <- rbind(circle1,circle[l-1,],circle[l,])}
    }

    circle1$x <- c2x(circle1$lon,circle1$lat)
    circle1$y <- c2y(circle1$lat)
    filling <- which(is.na(circle1))
    if(length(filling!=0)){
      lines(circle1$x,circle1$y, lwd=0.8)
    }
    else if(fill_A==FALSE) {lines(circle1$x,circle1$y, lwd=0.8)}
    else {polygon(circle1$x,circle1$y, lwd=0.8,col=rgb(1,0.9,0,0.25))}
  }

  #plot map
  if(on_plot==FALSE) {
    Map_KVII(grid=grid,center=center)
  }else{cat("
Double check the center meridian of the map!
")}


  #plot points from file
  if(S_file==TRUE){
    cat("Select file in .csv format")
    dat <- read.csv(file.choose())
    for(i in 1:nrow(dat)){
      S_lon <- dat[i,1]-center
      S_lon <- ifelse(S_lon<0,S_lon+360,S_lon)
      S_lon <- ifelse(S_lon>180,S_lon-360,S_lon)
      S_lat <- dat[i,2]
      if(A95==TRUE){
        circ <- dat[i,3]
        circle(S_lon,S_lat,circ,fill_A)
      }
      #select symbol
      if(symbol=="c") pch <- 21
      if(symbol=="s") pch <- 22
      if(symbol=="d") pch <- 23
      if(symbol=="t") pch <- 24

      site_x <- c2x(S_lon,S_lat)
      site_y <- c2y(S_lat)
      points(x=site_x,y=site_y,pch=pch, col="black",bg=col)
    }
  }
  if(S_file==FALSE){
    repeat{
      #plot point
      S_lon <- as.numeric(readline("Longitude -or press enter to exit-: "))
      if(is.na(S_lon)==TRUE) break
      S_lon <- S_lon-center
      S_lon <- ifelse(S_lon>180,
                      S_lon-360,S_lon)
      S_lon <- ifelse(S_lon<(-180),S_lon+360,S_lon)
      S_lat <- as.numeric(readline("Latitude: "))
      circ <- as.numeric(readline("Semi-angle of circle -press enter if no confidence angle is required-: "))
      if(is.na(circ)==FALSE){circle(S_lon,S_lat,circ,fill_A)}
      symbol <- readline("Symbol -press enter for circle-: ")
      #select symbol
      if(symbol=="") sym <- 21
      if(symbol=="c") sym <- 21
      if(symbol=="s") sym <- 22
      if(symbol=="d") sym <- 23
      if(symbol=="t") sym <- 24

      site_x <- c2x(S_lon,S_lat)
      site_y <- c2y(S_lat)
      color <- readline("Color -press enter for red-:")
      if(color=="") color <- "red"
      points(x=site_x,y=site_y,pch=sym, col="black",bg=color)
    }
  }
  if(export==TRUE){
    save_pdf(name="Map.pdf")
    cat("
Figure saved as Map.pdf
")
  }
}

#function that generate and plot confidence ellipses calculated from A95 back to directions, delta Inc > delta dec
generate_ellips <- function(D,I,delta_dec,delta_inc,col_d="red",col_u="white",col_l="black", symbol="c", on_plot=FALSE, save=FALSE, name="confidence_ellipse"){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #create circle with dummy bedding only for strain_DI to work
  circle <- as.data.frame(seq(0,360,2))
  circle$inc <- rep(90-delta_inc)
  colnames(circle) <- c("circ_dec","circ_inc")
  circle$dummy_az <- rep(0)
  circle$dummy_pl <- rep(0)
  #calculate parameter for adjusting delta_dec
  fol <- tan(d2r(delta_inc))/tan(d2r(delta_dec))
  #create matrix deforming circle
  M <- matrix_maker(Fol = fol,v1d = D,v1i = 0,v2d = 0,v2i = 90,v3d = ((D-90)%%360),v3i = 0, return_P=F)
  ell <- strain_DI(DIAP = circle,M = M)
  ell <- ell[,1:2]
  #uses bedding correction to place final ellipses in the right position
  ellipses <- bed_DI(ell,in_file = FALSE,bed_az = D,bed_plunge = 90-I,export = FALSE)
  colnames(ellipses) <- c("dec","inc")
  ellipses$x <- a2cx(abs(ellipses$inc),ellipses$dec)
  ellipses$y <- a2cy(abs(ellipses$inc),ellipses$dec)
  #restore screen
  par(fig=c(0,1,0,1))
  #standalone graph or on existing graph
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    equalarea()
  }
  #check for negative inclination
  inc <- I
  dec <- D
  UD <- ifelse(inc>0,"D","U")
  inc <- abs(inc)
  X <- a2cx(inc,dec)
  Y <- a2cy(inc,dec)
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}
  if(UD=="D"){
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg= col_d)
  }else{
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg=col_u)
  }
  lines(ellipses$x,ellipses$y,lty=1, col=col_l, lwd=1.8)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function that take data Dec_Inc and return the average Inc, E, and E declination
inc_E_finder <- function(DI, export=FALSE, name="I_E_Edec") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)

  #force V1 to positive
  V1inc <- abs(V1inc)
  V1dec <- ifelse(V1inc<0,(V1dec+180)%%360,V1dec)

  V2dec <- r2d(atan2(T_vec[2,2],T_vec[1,2]))
  V2dec <- ifelse(V2dec<0,V2dec+360,V2dec)
  V2inc <- r2d(asin(T_vec[3,2]/(sqrt((T_vec[1,2]^2)+(T_vec[2,2]^2)+(T_vec[3,2]^2)))))

  #Calculate difference between V1 and V2 to have the declination of V2 with respect to V1
  DV1V2 <- V1dec-V2dec
  DV1V2 <- ifelse(DV1V2<0,DV1V2+360,DV1V2)
  DV1V2 <- ifelse(DV1V2>90,ifelse(DV1V2<270,DV1V2-180,DV1V2),DV1V2)
  DV1V2 <- ifelse(DV1V2>270,DV1V2-360,DV1V2)

  E <- T_val[2]/T_val[3]

  inc_E <- as.data.frame(cbind(V1inc,E,DV1V2))
  #export result if requested
  if(export==TRUE){write.csv(round(inc_E,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(inc_E)
}

#Arason and Levi(2010) inclination only calculation
#adepted from the original fortran source code ARALEV available at http://hergilsey.is/arason/paleomag/aralev.txt
inc_only <- function(DI,dec=TRUE, print=TRUE,export=TRUE, name="Inclination_only",return=TRUE, arith_stat=FALSE) {
  #The Arason-Levi MLE Iteration Formula 1
  AL1 <- function(th, n, the, ak) {
    dr <- 0.0174532925199433  # Degrees to radians (pi/180)
    s <- 0
    c <- 0

    for (i in 1:n) {
      x <- ak * sin(the * dr) * sin(th[i] * dr)
      bessel_result <- bessel(x)
      bi1i0 <- bessel_result[3]

      s <- s + sin(th[i] * dr) * bi1i0
      c <- c + cos(th[i] * dr)
    }

    AL1 <- atan2(s, c) / dr
    AL1 <- ifelse(AL1 < 0.000001, 0.000001, AL1)
    AL1 <- ifelse(AL1 > 179.999999, 179.999999, AL1)

    return(AL1)
  }
  #The Arason-Levi MLE Iteration Formula 2
  AL2 <- function(th, n, the, ak) {
    dr <- 0.0174532925199433  # Degrees to radians (pi/180)
    dn <- n

    s <- 0
    c <- 0
    for (i in 1:n) {
      x <- ak * sin(the * dr) * sin(th[i] * dr)
      Btemp <- bessel(x)
      bi1i0 <- Btemp[3]

      s <- s + sin(th[i] * dr) * bi1i0
      c <- c + cos(th[i] * dr)
    }

    x <- dn * (1 / tanh(ak)) - cos(the * dr) * c - sin(the * dr) * s
    AL2 <- 1e10
    if (x / dn > 1e-10) {
      AL2 <- dn / x
    }
    if (AL2 < 1e-06) {
      AL2 <- 1e-06
    }

    return(AL2)
  }
  #Evaluation of the Hyperbolic Bessel functions I0(x), I1(x) and their ratio I1(x)/I0(x).
  bessel <- function(x) {
    p <- c(1.0, 3.5156229, 3.0899424, 1.2067492, 0.2659732, 0.360768e-1, 0.45813e-2)
    q <- c(0.39894228, 0.1328592e-1, 0.225319e-2, -0.157565e-2, 0.916281e-2, -0.2057706e-1, 0.2635537e-1, -0.1647633e-1, 0.392377e-2)
    u <- c(0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1, 0.301532e-2, 0.32411e-3)
    v <- c(0.39894228, -0.3988024e-1, -0.362018e-2, 0.163801e-2, -0.1031555e-1, 0.2282967e-1, -0.2895312e-1, 0.1787654e-1, -0.420059e-2)

    if (abs(x) < 3.75) {
      t <- (x / 3.75)^2
      b0 <- p[1] + t * (p[2] + t * (p[3] + t * (p[4] + t * (p[5] + t * (p[6] + t * p[7])))))
      b1 <- x * (u[1] + t * (u[2] + t * (u[3] + t * (u[4] + t * (u[5] + t * (u[6] + t * u[7]))))))
      bi0e <- b0 / exp(abs(x))
      bi1e <- b1 / exp(abs(x))
      bi1i0 <- b1 / b0
    } else {
      t <- 3.75 / abs(x)
      b0 <- q[1] + t * (q[2] + t * (q[3] + t * (q[4] + t * (q[5] + t * (q[6] + t * (q[7] + t * (q[8] + t * q[9])))))))
      b1 <- v[1] + t * (v[2] + t * (v[3] + t * (v[4] + t * (v[5] + t * (v[6] + t * (v[7] + t * (v[8] + t * v[9])))))))
      if (x < 0) b1 <- -b1
      bi0e <- b0 / sqrt(abs(x))
      bi1e <- b1 / sqrt(abs(x))
      bi1i0 <- b1 / b0
    }

    return(c(bi0e, bi1e, bi1i0))
  }
  #Evaluation of the Hyperbolic Cotangens function coth(x)
  coth <- function(x) {
    if (x == 0) {
      return(0)
    }

    t <- abs(x)

    if (t < 0.001) {
      return(1 / t + t / 3 - t^3 / 45 + 2 * t^5 / 945)
    } else if (t <= 15) {
      ep <- exp(t)
      em <- exp(-t)
      return((ep + em) / (ep - em))
    } else {
      return(1)
    }

    if (x < 0) {
      return(-coth)
    }
  }
  #Evaluation of the Log-Likelihood function for inclination-only data.
  xlik <- function(th, n, the, ak) {
    dr <- 0.0174532925199433         # Degrees to radians (pi/180)
    pi <- 180 * dr
    dn <- n

    # Illegal use
    if (n < 1) {
      cat("ERROR: Data missing in xlik\n")
      return(-1e10)
    }

    if (n > 10000) {
      cat("ERROR: Too small dimension in xlik\n")
      return(-1e10)
    }

    # Uncomment the following lines if you want to check the range of ak
    # if (ak < 0) {
    #   cat("ERROR: Out of range in xlik\n")
    #   return(-1e10)
    # }

    # A1(k) = N ln(k) - N ln(sinh k) - N ln(2)
    a1 <- 0

    if (ak >= 0 && ak < 0.01) {
      q <- -ak * (1 - ak * (2 / 3 - ak * (1 / 3 - ak * (2 / 15 - ak * (8 / 45)))))
      a1 <- dn * (-log(2) - log(1 + q) - ak)
    } else if (ak >= 0.01 && ak <= 15) {
      a1 <- dn * (log(ak) - log(1 - exp(-2 * ak)) - ak)
    } else {
      a1 <- dn * (log(ak) - ak)
    }

    # A2(k,t,ti) = Sum(k cos t cos ti) + Sum(ln(BessIo(k sin t sin ti)))
    a2 <- 0

    for (i in 1:n) {
      x <- ak * sin(the * dr) * sin(th[i] * dr)
      bi0e <- bessel(x)[[1]]
      bi1e <- bessel(x)[[2]]
      bi1i0 <- bessel(x)[[3]]
      a2 <- a2 + ak * cos((th[i] - the) * dr) + log(bi0e)
    }

    # A3(ti) = Sum( ln(sin(ti)) ), Note: 0.000001 < ti < 179.999999
    a3 <- 0

    for (i in 1:n) {
      x <- th[i]
      if (x < 1e-6) x <- 1e-6
      if (x > 179.999999) x <- 179.999999
      a3 <- a3 + log(sin(x * dr))
    }

    # The log-likelihood function
    xlik <- a1 + a2 + a3

    return(xlik)
  }
  #Calculation of the arithmetic mean of inclination-only data
  armean <- function(xinc) {
    dr <- 0.01745329252         # Degrees to radians (pi/180)
    t63max <- 105.070062145     # 63 % of a sphere.
    a95max <- 154.158067237     # 95 % of a sphere.
    dn <- length(xinc)

    s <- sum(xinc)
    s2 <- sum(xinc^2)

    ainc <- s / dn

    sd <- 0
    ak <- -1

    if (dn > 1) {
      sd <- sqrt((s2 - s^2 / dn) / (dn - 1))
      ak <- (dn - 1) / ((s2 - s^2 / dn) * dr^2)
    }

    nf <- dn - 1

    tval_63 <- qt(0.63, df = nf)
    t63 <- tval_63 * sd

    tval_95 <- qt(0.95, df = nf)
    a95 <- tval_95 * sd / sqrt(dn)

    result <- as.data.frame(matrix(ncol=5, nrow=1))
    result[1] <- dn
    result[2] <- round(ainc, digits=2)
    result[3] <- round(ak, digits=2)
    result[4] <- round(t63,digits = 2)
    result[5] <- round(a95, digits = 2)
    colnames(result) <- c("N","Inc","Precision","Angular st.dev(63%)","a95")

    return(result)
  }

  #isolate inclination if file contains also declination
  DI <- na.omit(DI)
  if(dec==TRUE) {DI <- DI[,1:2]
  }else if(dec==F) {DI <- DI[,1]}
  if(dec==TRUE) {
    colnames(DI) <- c("dec","inc")
    xinc <- DI$inc
  }else{
    #convert the imported file in a list of numbers
    inc <- as.list(inc)
    xinc <- inc[[1]]
  }
  #generate only arithmetic statistic if requested
  if(arith_stat==TRUE){
    result <- armean(xinc)
    #print if request
    if(export==TRUE){write.csv(result,paste(name,".csv"), row.names = F)}

    if(print==TRUE){
      cat("Arithmetic average inclination result:

N:", result[1,1],"
Inclination:", result[1,2],"
Precision:",result[1,3],"
St.Dev_63:", result[1,4],"
alpha_95:", result[1,5],"

")
    }
    if(return==TRUE) {return(result)}
  }else{
    # main routine
    # Degrees to radians (pi/180)
    dr <- 0.0174532925199433
    # 63 % of a sphere
    t63max <- 105.070062145
    # 95 % of a sphere
    a95max <- 154.158067237
    n <- length(xinc)
    th <- numeric(n)
    dn <- n
    ierr <- 1
    #sdata checking and warnings
    if (n == 1) {
      stop("Only one observed inclination\n", call.=F)
    }
    if (n > 10000) {
      assign("inc_warn","Too many directions", envir = .GlobalEnv)
      stop("Too many directions, max=10000\n", call.=F)
    }
    if(length(unique(xinc))==1){
      assign("inc_warn","Directions are all identical", envir = .GlobalEnv)
      stop("Directions are all identical\n", call.=F)
    }
    if (any(xinc>90) | any(xinc<(-90))) {
      assign("inc_warn","Inclination must be between -90 and 90", envir = .GlobalEnv)
      stop("Inclination must be between -90 and 90\n", call.=F)
    }
    #file with co-incl
    for (i in 1:n) {
      th[i] <- 90 - xinc[i]
    }

    s <- sum(th)
    s2 <- sum(th^2)
    c <- sum(cos(th * dr)) / dn
    # initial theta guess (17)
    rt <- s / dn
    x <- (s2 - s^2 / dn) * dr^2
    rk <- ifelse(x / (dn - 1) > 1e-10, (dn - 1) / x, 1e10)
    rt1 <- rt
    rk1 <- rk
    ie1 <- 0

    the1 <- rt
    akap1 <- rk
    for (j in 1:10000) {
      rt <- AL1(th, n, rt, rk)
      rk <- AL2(th, n, rt, rk)
      dt <- abs(rt - the1)
      dk <- abs((rk - akap1) / rk)
      if (j > 10 && dt < 1e-6 && dk < 1e-6) break
      the1 <- rt
      akap1 <- rk
    }
    #likelihood for theta e k
    #ie1 <- 0
    the1 <- rt
    akap1 <- rk
    xl1 <- xlik(th, n, rt, rk)

    #likelihood for theta=0
    rt <- 0
    rk <- rk1
    akap2 <- rk
    ie2 <- 0
    for (j in 1:10000) {
      x <- coth(rk) - c
      if (x > 1e-10) {
        rk <- 1 / x
      } else {
        rk <- 1e10
      }
      dk <- abs((rk - akap2) / rk)
      if (j > 4 && dk < 1e-6) break
      if (rk < 1e-6) break
      akap2 <- rk
    }
    ie2 <- 1
    the2 <- 0
    akap2 <- rk
    xl2 <- xlik(th, n, rt, rk)

    #likelihood for theta=180
    rt <- 180
    rk <- rk1
    akap3 <- rk
    ie3 <- 0
    for (j in 1:10000) {
      x <- coth(rk) + c
      if (x > 1e-10) {
        rk <- 1 / x
      } else {
        rk <- 1e10
      }
      dk <- abs((rk - akap3) / rk)
      if (j > 4 && dk < 1e-6) break
      if (rk < 1e-6) break
      akap3 <- rk
    }
    ie3 <- 1
    the3 <- 180
    akap3 <- rk
    xl3 <- xlik(th, n, rt, rk)

    #likelihood for k=0
    rt <- 90
    rk <- 0
    the4 <- rt
    akap4 <- rk
    xl4 <- xlik(th, n, rt, rk)

    isol <- 1
    ierr <- ie1
    if (xl2 > xl1) {
      the1 <- the2
      akap1 <- akap2
      xl1 <- xl2
      isol <- 2
      ierr <- 1
    }

    #compares solutions
    if (xl3 > xl1) {
      the1 <- the3
      akap1 <- akap3
      xl1 <- xl3
      isol <- 3
      ierr <- 1
    }

    if (xl4 > xl1) {
      the1 <- the4
      akap1 <- akap4
      xl1 <- xl4
      isol <- 4
      ierr <- 0
    }

    ainc <- 90 - the1
    ak <- akap1
    if (ierr != 0) {
      assign("inc_warn","Convergence problems", envir = .GlobalEnv)
      cat("Convergence problems\n")
    }

    #Test of robustness with 16 surrounding points
    for (i in 1:16) {
      x <- i
      rt <- the1 + 0.01 * cos(22.5 * x * dr)
      if (rt < 0 || rt > 180) break
      rk <- akap1 * (1 + 0.001 * sin(22.5 * x * dr))
      xl <- xlik(th, n, rt, rk)
      if (xl > xl1) {
        ierr <- ierr + 2
        assign("inc_warn","Robustness failure", envir = .GlobalEnv)
        cat("Robustness failure\n")
      }
    }

    if (akap1 >= 20) {
      co <- 1 + log(1 - 0.63) / akap1
    } else if (akap1 > 0.1 && akap1 < 20) {
      co <- 1 + log(1 - 0.63 * (1 - exp(-2 * akap1))) / akap1
    } else if (akap1 <= 0.1) {
      co <- -0.26 + 0.4662 * akap1
    }

    #calculate confidences
    t63 <- 90 - (90*sign(co))
    if (abs(co) < 1) {
      t63 <- 90 - atan(co / sqrt(1 - co^2)) / dr
    }
    if (t63 > t63max) {
      t63 <- t63max
    }

    co <- 1 - (dn - 1) * (20^(1 / (dn - 1)) - 1) / (dn * (akap1 - 1) + 1)
    a95 <- 90 - (90*sign(co))
    if (abs(co) < 1) {
      a95 <- 90 - atan(co / sqrt(1 - co^2)) / dr
    }
    if (a95 > a95max) {
      a95 <- a95max
    }

    #calculates arith mean
    ari_mean <- armean(xinc)

    #compile result file
    result <- as.data.frame(matrix(ncol=5, nrow=1))
    result[1] <- n
    result[2] <- round(ainc, digits=2)
    result[3] <- round(ak, digits=2)
    result[4] <- round(t63,digits = 2)
    result[5] <- round(a95, digits = 2)
    result[6] <- round(ari_mean[1,2],digits = 2)
    colnames(result) <- c("N","Inc","Precision","Angular st.dev(63%)","a95","Aritm. mean")

    #print if request
    if(export==TRUE){write.csv(result,paste(name,".csv"), row.names = F)}

    if(print==TRUE){
      cat("Arason-Levi inclination only result:

N:", result[1,1],"
Inclination:", result[1,2],"
Precision:",result[1,3],"
St.Dev_63:", result[1,4],"
alpha_95:", result[1,5],"
Aritm. mean:",result[1,6],"

")
    }
  }
  if(return==TRUE) {return(result)}
}

#plot equal area of Arason and Levi(2010) inclination only calculation
inc_plot <- function(DI,dec=TRUE,plot=TRUE,bimodal=FALSE,on_plot=TRUE, col="black", print=TRUE,export=TRUE, save=TRUE,name="Inc_only", arith_stat=FALSE,Shiny=FALSE){
  #import dplyr for filter_ALL
  library(dplyr)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  DI <- na.omit(DI)
  #splits modes
  if(Shiny==TRUE){
    if(arith_stat==FALSE){
      inconly_stat <- data.frame(matrix(nrow=0,ncol=6))
      colnames(inconly_stat) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
    }else{
      inconly_stat <- data.frame(matrix(nrow=0,ncol=5))
      colnames(inconly_stat) <- c("N","Inc","x","St.D.","a95")
    }
  }
  if(bimodal==TRUE){
     if(dec==TRUE) {DI <- DI[,1:2]
    }else if(dec==F) {DI <- DI[,1]}
    dirs <- DI
    ifelse(dec==TRUE, colnames(dirs) <- c("dec","inc"), colnames(dirs) <- "inc")
    dirs_D <- filter_all(dirs, all_vars(inc>0))
    dirs_U <- filter_all(dirs, all_vars(inc<=0))
    #down_pointing
    if(print==TRUE){cat("Down-pointing\n")}
    inc_stat_D <- inc_only(DI = dirs_D,dec = dec, print = print,export=export, name=paste(name,"_down"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_D) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_D) <- "Mode 1"
        inconly_stat <- rbind(inconly_stat,inc_stat_D)
      }else{
        colnames(inc_stat_D) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_D) <- "Mode 1"
        inconly_stat <- rbind(inconly_stat,inc_stat_D)
      }
    }
    #up_pointing
    if(print==TRUE){cat("Up-pointing\n")}
    inc_stat_U <- inc_only(DI = dirs_U,dec = dec, print = print,export=export, name=paste(name,"_up"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_U) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_U) <- "Mode 2"
        inconly_stat <- rbind(inconly_stat,inc_stat_U)
      }else{
        colnames(inc_stat_U) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_U) <- "Mode 2"
        inconly_stat <- rbind(inconly_stat,inc_stat_U)
      }
    }
    #all_down_pointing
    if(print==TRUE){cat("All data\n")}
    dirs$inc <- abs(dirs$inc)
    inc_stat_ALL <- inc_only(DI = dirs,dec = dec, print = print,export=export, name=paste(name,"_all"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }else{
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }
    }
  }else{
    dirs <- DI
    inc_stat_ALL <- inc_only(DI = dirs,dec = dec, print = print,export=export, name=paste(name,"_all"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }else{
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }
    }
  }
  #plot if requested
  if(plot==TRUE){
    if(bimodal==TRUE){
      #create circles
      circle_D <- as.data.frame(seq(0,90,1))
      colnames(circle_D) <- "dec"
      circle_D$inc <- rep(inc_stat_D$Inc)
      #create a column with horizontal inclination for doubling the confidence area if crosses 0 inc
      if(inc_stat_D$Inc-inc_stat_D$a95 < 0){
        circle_D$x_hrz <- a2cx(0,circle_D$dec)
        circle_D$y_hrz <- a2cy(0,circle_D$dec)
      }
      circle_D$inc_l <- rep(abs(inc_stat_D$Inc-inc_stat_D$a95))
      circle_D$inc_h <- ifelse((inc_stat_D$Inc+inc_stat_D$a95)<90, rep(inc_stat_D$Inc+inc_stat_D$a95),rep(90))
      circle_D$x <- a2cx(circle_D$inc,circle_D$dec)
      circle_D$y <- a2cy(circle_D$inc,circle_D$dec)
      circle_D$x_l <- a2cx(circle_D$inc_l,circle_D$dec)
      circle_D$y_l <- a2cy(circle_D$inc_l,circle_D$dec)
      circle_D$x_h <- a2cx(circle_D$inc_h,circle_D$dec)
      circle_D$y_h <- a2cy(circle_D$inc_h,circle_D$dec)

      circle_U <- as.data.frame(seq(270,360,1))
      colnames(circle_U) <- "dec"
      circle_U$inc <- rep(abs(inc_stat_U$Inc))
      #create a column with horizontal inclination for doubling the confidence area if crosses 0 inc
      if(abs(inc_stat_U$Inc)-inc_stat_U$a95 < 0){
        circle_U$x_hrz <- a2cx(0,circle_U$dec)
        circle_U$y_hrz <- a2cy(0,circle_U$dec)
      }
      circle_U$inc_l <- rep(abs(abs(inc_stat_U$Inc)-inc_stat_U$a95))
      circle_U$inc_h <- ifelse((abs(inc_stat_U$Inc)+inc_stat_U$a95)<90, rep(abs(inc_stat_U$Inc)+inc_stat_U$a95), rep(90))
      circle_U$x <- a2cx(circle_U$inc,circle_U$dec)
      circle_U$y <- a2cy(circle_U$inc,circle_U$dec)
      circle_U$x_l <- a2cx(circle_U$inc_l,circle_U$dec)
      circle_U$y_l <- a2cy(circle_U$inc_l,circle_U$dec)
      circle_U$x_h <- a2cx(circle_U$inc_h,circle_U$dec)
      circle_U$y_h <- a2cy(circle_U$inc_h,circle_U$dec)
      circle_ALL <- as.data.frame(seq(90,270,1))
      colnames(circle_ALL) <- "dec"
    }else{
      #if biomdal is false draw a complete ALL circle
      circle_ALL <- as.data.frame(seq(0,360,1))
      colnames(circle_ALL) <- "dec"
    }
    #complete circle with all even if not bimodal
    circle_ALL$inc <- rep(abs(inc_stat_ALL$Inc))
    #create a column with horizontal inclination for doubling the confidence area if crosses 0 inc
    if(abs(inc_stat_ALL$Inc)-inc_stat_ALL$a95 < 0){
      circle_ALL$x_hrz <- a2cx(0,circle_ALL$dec)
      circle_ALL$y_hrz <- a2cy(0,circle_ALL$dec)
    }
    circle_ALL$inc_l <- rep(abs(abs(inc_stat_ALL$Inc)-inc_stat_ALL$a95))
    circle_ALL$inc_h <- ifelse((abs(inc_stat_ALL$Inc)+inc_stat_ALL$a95)<90,rep(abs(inc_stat_ALL$Inc)+inc_stat_ALL$a95),rep(90))
    circle_ALL$x <- a2cx(circle_ALL$inc,circle_ALL$dec)
    circle_ALL$y <- a2cy(circle_ALL$inc,circle_ALL$dec)
    circle_ALL$x_l <- a2cx(circle_ALL$inc_l,circle_ALL$dec)
    circle_ALL$y_l <- a2cy(circle_ALL$inc_l,circle_ALL$dec)
    circle_ALL$x_h <- a2cx(circle_ALL$inc_h,circle_ALL$dec)
    circle_ALL$y_h <- a2cy(circle_ALL$inc_h,circle_ALL$dec)

    #plot_circles
    if(on_plot==FALSE) equalarea()
    if(bimodal==TRUE){
      #Down
      lines(circle_D$x_l,circle_D$y_l,lty=2)
      lines(circle_D$x_h,circle_D$y_h,lty=2)
      #splits the confidence area in two if crosses zero inclination
      if(inc_stat_D$Inc-inc_stat_D$a95 < 0){
        conf1_D <- data.frame(cbind(circle_D$x_l,circle_D$y_l))
        confh_D <- data.frame(cbind(circle_D$x_hrz,circle_D$y_hrz))
        confh_D <- confh_D[nrow(confh_D):1,]
        first_D_area <- rbind(conf1_D,confh_D)
        polygon(first_D_area, col=rgb(0,0,1,0.30),border=NA)
        conf2_D <- data.frame(cbind(circle_D$x_h, circle_D$y_h))
        second_D_area <- rbind(conf2_D,confh_D)
        polygon(second_D_area, col=rgb(0,0,1,0.30),border=NA)

      }else{
        conf1_D <- data.frame(cbind(circle_D$x_l,circle_D$y_l))
        conf2_D <- data.frame(cbind(circle_D$x_h, circle_D$y_h))
        conf2_D <- conf2_D[nrow(conf2_D):1,]
        conf_D <- rbind(conf1_D,conf2_D)
        polygon(conf_D, col=rgb(0,0,1,0.30),border=NA)
      }
      lines(circle_D$x,circle_D$y,col=col, lwd=1.5)

      #Up
      lines(circle_U$x_l,circle_U$y_l,lty=2)
      lines(circle_U$x_h,circle_U$y_h,lty=2)
      #splits the confidence area in two if crosses zero inclination
      if(abs(inc_stat_U$Inc)-inc_stat_U$a95 < 0){
        conf1_U <- data.frame(cbind(circle_U$x_l,circle_U$y_l))
        confh_U <- data.frame(cbind(circle_U$x_hrz,circle_U$y_hrz))
        confh_U <- confh_U[nrow(confh_U):1,]
        first_U_area <- rbind(conf1_U,confh_U)
        polygon(first_U_area, col=rgb(0,0,1,0.30),border=NA)
        conf2_U <- data.frame(cbind(circle_U$x_h, circle_U$y_h))
        second_U_area <- rbind(conf2_U,confh_U)
        polygon(second_U_area, col=rgb(0,0,1,0.30),border=NA)
      }else{
        conf1_U <- data.frame(cbind(circle_U$x_l,circle_U$y_l))
        conf2_U <- data.frame(cbind(circle_U$x_h, circle_U$y_h))
        conf2_U <- conf2_U[nrow(conf2_U):1,]
        conf_U <- rbind(conf1_U,conf2_U)
        polygon(conf_U, col=rgb(0,1,1,0.30),border=NA)
      }
      lines(circle_U$x,circle_U$y,col=col, lwd=1.5 )
    }
    #ALL
    lines(circle_ALL$x_l,circle_ALL$y_l,lty=2)
    lines(circle_ALL$x_h,circle_ALL$y_h,lty=2)
    if(abs(inc_stat_ALL$Inc)-inc_stat_ALL$a95 < 0){
      conf1_ALL <- data.frame(cbind(circle_ALL$x_l,circle_ALL$y_l))
      confh_ALL <- data.frame(cbind(circle_ALL$x_hrz,circle_ALL$y_hrz))
      confh_ALL <- confh_ALL[nrow(confh_ALL):1,]
      first_ALL_area <- rbind(conf1_ALL,confh_ALL)
      polygon(first_ALL_area, col=rgb(0,0,1,0.30),border=NA)
      conf2_ALL <- data.frame(cbind(circle_ALL$x_h, circle_ALL$y_h))
      second_ALL_area <- rbind(conf2_ALL,confh_ALL)
      polygon(second_ALL_area, col=rgb(0,0,1,0.30),border=NA)

    }else{
      conf1_ALL <- data.frame(cbind(circle_ALL$x_l,circle_ALL$y_l))
      conf2_ALL <- data.frame(cbind(circle_ALL$x_h, circle_ALL$y_h))
      conf2_ALL <- conf2_ALL[nrow(conf2_ALL):1,]
      conf_ALL <- rbind(conf1_ALL,conf2_ALL)
      infill <- ifelse(inc_stat_ALL$Inc<0,rgb(0,1,1,0.30),rgb(0,0,1,0.30))
      if(bimodal==TRUE){infill <- rgb(1,0,0,0.30)}
      polygon(conf_ALL, col=infill,border=NA)
    }
    lines(circle_ALL$x,circle_ALL$y,col=col, lwd=1.5)
    if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
  }
  if(Shiny==TRUE){return(inconly_stat[,1:5])} #cut last column with Arithmetic mean
}

#create matrix from fol,lin, and dec inc of vectors
matrix_maker <- function(Fol=1,Lin=1,v1d,v1i,v2d,v2i,v3d,v3i, export=FALSE, name="matrix",return_P=TRUE){
  library(matlib)
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  P <- Fol*Lin
  #eigenvalues
  v1 <- (3*Lin)/(Lin+(1/Fol)+1)
  v2 <- v1/Lin
  v3 <- v1/P
  A <- as.data.frame(matrix(c(v1d,v1i,
                              v2d,v2i,
                              v3d,v3i),
                            ncol=2,byrow=TRUE))
  colnames(A) <- c("dec","inc")
  A$x <- cos(d2r(A$dec))*cos(d2r(A$inc))
  A$y <- sin(d2r(A$dec))*cos(d2r(A$inc))
  A$z <- sin(d2r(A$inc))
  A_vec <- t(A[,3:5])
  A_val <- c(v1,0,0,0,v2,0,0,0,v3)
  A_val <- matrix(A_val,nrow=3,byrow=TRUE)
  M <- A_vec%*%A_val%*%inv(A_vec)
  if(return_P==TRUE){
    cat(paste("P:",P,"
"))
  }
  #export matrix if requested
  if(export==TRUE){write.csv(round(M,digits=5),paste(name,".csv"),row.names = TRUE)}
  return(round(M, digits=5))
}

#function that plots KavrayskiyVII geographic projection
Map_KVII <- function(grid=30, center=0, title="",seaCol="light cyan",landCol="light green",gridCol="gray") {
  library(rlist)
  if(center>180 | center<(-180)) stop("Please set center between -180° and 180°",call. = F)
  if(grid>90) stop("Please set the grid between 1° and 90°. If grid=0, grid is not plotted",call. = F)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {((3*d2r(lon))/2)*(sqrt((1/3)-((d2r(lat)/pi)^2)))}
  c2y <- function(lat) {d2r(lat)}

  #fix frame
  plot(NA, xlim=c(-2.6,2.6), ylim=c(-1.5,1.5), asp=1,
       xlab=title, xaxt="n",ylab="", yaxt="n", axes=FALSE)
  #import coastline from PmagDiR
  if (abs(center)==180){
    cl <- world_coastline_180
  }else{cl <- world_coastline}
  # set the coastline offset if longitude is not 0
  if(center!=0 | abs(center)!=180){
    #create the coastline-breaks file (breaks separating the single continental lines)
    sep <- as.data.frame(1)
    colnames(sep) <- "index"
    #finds all breaks in the coastline file and index them
    list <- as.data.frame(which(is.na(cl[,1]),arr.ind = TRUE))
    colnames(list) <- "index"
    sep <- rbind(sep, list)
    rm(list)
    #create empty list of all continents coastline
    conts <- list()
    #function that isolate coastlines and append to list
    for(i in 2:nrow(sep)){
      contour <- cl[sep[i-1,1]:sep[i,1],]
      contour <- na.omit(contour)
      conts <- list.append(conts,contour)
    }
    new_cl <- as.data.frame(matrix(nrow = 0,ncol = 2))
    colnames(new_cl) <- c("lon","lat")
    #for every coastline it changes the longitude depending on the new center and it fix it (-180<l<+180)
    for(i in 1:length(conts)){
      clc <- conts[[i]]
      clc$lon <- clc$lon-center
      clc$lon <- ifelse(clc$lon<(-180),clc$lon+360,clc$lon)
      clc$lon <- ifelse(clc$lon>180,clc$lon-360,clc$lon)
      #creates a new dataframe for single coastlines with breaks
      clc1 <- as.data.frame(matrix(nrow = 0,ncol = 2))
      colnames(clc1) <- c("lon","lat")
      #creates the breaking line
      breaker <- as.data.frame(t(c(NA,NA)))
      colnames(breaker) <- c("lon","lat")
      #when two points are on the different side of the map, based on 350° distance, it put break in between
      for(l in 2:nrow(clc)){
        if(abs(clc[l-1,1]-clc[l,1])>350){
          provv1 <- clc[l-1,]
          provv2 <- clc[l,]
          clc1 <- rbind(clc1,provv1,breaker,provv2)
        }else{clc1 <- rbind(clc1,clc[l-1,],clc[l,])}
      }
      #puts break after the new continent line
      clc1 <- rbind(clc1,breaker)
      #appends all new continent lines
      new_cl <- rbind(new_cl,clc1)
    }
  }
  #plot sea
  bord_left <- as.data.frame(matrix(ncol = 2, nrow = 181))
  bord_left[,2] <- -90:90
  bord_left[,1] <- rep(-180)
  bord_right <- as.data.frame(matrix(ncol = 2, nrow = 181))
  bord_right[,2] <- 90:-90
  bord_right[,1] <- rep(180)
  bord <- rbind(bord_left,bord_right)
  bord$x <- c2x(bord[,1], bord[,2])
  bord$y <- c2y(bord[,2])
  if(center==0 | abs(center)==180){
    polygon(bord$x,bord$y, col=seaCol,border = NA)
  }

  #set coastline if longitude is greenwich
  if(center==0 | abs(center)==180){
    new_cl <- cl
    polygon(x = c2x(new_cl$lon,new_cl$lat),
                        y = c2y(new_cl$lat), col=landCol, border=landCol)
    }else{
      lines(x = c2x(new_cl$lon,new_cl$lat),
        y = c2y(new_cl$lat), col="black")
    }

  #plot grid only if different from 0
  if(grid!=0){
    #plot_main_parallel
    #longitude circle
    lats <- seq(-(90-grid),(90-grid),grid)
    for(i in lats){
      lon_lat_p <-  as.data.frame(-180:180)
      lon_lat_p$lat <- rep(i)
      lon_lat_p$x <- c2x(lon_lat_p[,1],lon_lat_p[,2])
      lon_lat_p$y <- c2y(lon_lat_p[,2])
      lines(lon_lat_p$x,lon_lat_p$y,col=gridCol, pch=16, cex=0.3, lty=1)
    }
    #plot_main_meridians
    #fix meridians if center is not greenwich
    Gr <- (-center)
    LonLeft <- seq((Gr),-180,-grid)
    LonRight <- seq(Gr,180,grid)
    LonRight <- LonRight[-1]
    #plot meridians left of Greenwich
    for(i in LonLeft){
      lat_lon_m <- as.data.frame(seq(-89,89,1))
      lat_lon_m$lon <- rep(i)
      lat_lon_m$x <- c2x(lat_lon_m[,2],lat_lon_m[,1])
      lat_lon_m$y <- c2y(lat_lon_m[,1])
      lines(lat_lon_m$x,lat_lon_m$y,col=gridCol, pch=16, cex=0.3, lty=1)
    }
    #plot meridians right of Greenwich
    for(i in LonRight){
      lat_lon_m <- as.data.frame(seq(-89,89,1))
      lat_lon_m$lon <- rep(i)
      lat_lon_m$x <- c2x(lat_lon_m[,2],lat_lon_m[,1])
      lat_lon_m$y <- c2y(lat_lon_m[,1])
      lines(lat_lon_m$x,lat_lon_m$y,col=gridCol, pch=16, cex=0.3, lty=1)
    }
  }
  #plot contour
  polygon(bord$x,bord$y, col=NA)
}

#calculate PCA-derived direction and MAD from demagnetization steps
PCA_DI <- function(DII,anchor="f", export=FALSE,name="PCA") {
  #degree to radians and VV
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DII
  colnames(data) <- c("dec", "inc","int")
  #directions in Cartesian coordinates
  data$x <- data$int*cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- data$int*sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- data$int*sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #copy coordinates for anchored directions
  if (anchor=="a"){
    data$xn <- data$x
    data$yn <- data$y
    data$zn <- data$z
  } else if(anchor=="i") {
    #includes origin and calculate new center of mass
    newrow <- c(0,0,0,0,0,0)
    data <- rbind(newrow,data)
    data$xn <- data$x-x_av
    data$yn <- data$y-y_av
    data$zn <- data$z-z_av
  } else {
    #calculate coordinates with new center of mass
    data$xn <- data$x-x_av
    data$yn <- data$y-y_av
    data$zn <- data$z-z_av
  }
  #elements of the distribution matrix
  T_elements <- c(sum((data$xn)*(data$xn)),sum(data$xn*data$yn),sum(data$xn*data$zn),
                  sum(data$yn*data$xn),sum(data$yn*data$yn),sum(data$yn*data$zn),
                  sum(data$zn*data$xn),sum(data$zn*data$yn),sum(data$zn*data$zn))

  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$value
  #calculate dec inc of max variance
  Vdec <- (r2d(atan2(T_vec[2,1],T_vec[1,1])))%%360
  Vinc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  MAD <- r2d(atan(sqrt(((T_val[2])+(T_val[3]))/T_val[1])))
  N <- length(data[,1])

  dirs <- cbind(Vdec,Vinc,MAD,N)
  colnames(dirs) <- c("Dec", "Inc","MAD","N")
  if(export==TRUE){write.csv(round(dirs,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(dirs)
}

#function that plots a95 with dec,inc,a95
plot_a95 <- function(D,I,a, col_d="red",col_u="white",col_l="black", symbol="c", on_plot=FALSE, save=FALSE, name="F_a95"){
  library("dplyr", warn.conflicts = FALSE)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #functions spherical (Dec=x, Inc=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  #save declination and inc and calculate new system for rotation
  dec <- D
  newSdec <- ifelse((D+180)>360,D-180,D+180)
  inc <- I
  newSinc <- 90-I
  newSdecr <- d2r(newSdec)
  newSincr <- d2r(newSinc)
  a95 <- a
  circle <- as.data.frame(matrix(ncol=2,nrow=0))
  #loop that create a95 and rotate it around new coordinate (dec, inc)
  for (i in seq(0,360,2)){
    circleP <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(i,(90-a95))
    y <- s2cy(i,(90-a95))
    z <- s2cz(90-a95)
    vec <- as.matrix(c(x,y,z))
    R_elements <- c(cos(newSincr)*cos(newSdecr), -sin(newSdecr), -sin(newSincr)*cos(newSdecr),
                    cos(newSincr)*sin(newSdecr), cos(newSdecr), -sin(newSincr)*sin(newSdecr),
                    sin(newSincr), 0, cos(newSincr))
    R <- matrix(R_elements,nrow=3, byrow=TRUE)
    newvec <- R%*%vec
    newdec <- r2d(atan2(newvec[2,1],newvec[1,1]))
    newdec <- ifelse(newdec<0,newdec+360,newdec)
    #absolute value avoid point outside the graph
    newinc <- abs(r2d(asin(newvec[3,1])))
    circleP[1,1:2] <- c(newdec,newinc)
    circle <- rbind(circle,circleP)
  }
  colnames(circle) <- c("dec","inc")
  circle$x <- a2cx(circle$inc,circle$dec)
  circle$y <- a2cy(circle$inc,circle$dec)
  #restore screen
  par(fig=c(0,1,0,1))
  #standalone graph or on existing graph
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    equalarea()
    }
  UD <- ifelse(inc>0,"D","U")
  inc <- abs(inc)
  X <- a2cx(inc,dec)
  Y <- a2cy(inc,dec)
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  if(UD=="D"){
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg= col_d)
  }else{
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg=col_u)
  }
  lines(circle$x,circle$y,lty=1, col=col_l, lwd=1.8)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot A95 on a spherical orthographic plot
plot_PA95 <- function(lon,lat,A,lon0=0,lat0=90,grid=30, col_f="red",col_b="white",col_l="black",col_A=rgb(1,0,0,0.30), symbol="c",size=1, coast=FALSE, on_plot=FALSE, save=FALSE, name="A95"){
  library("dplyr", warn.conflicts = FALSE)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}

  #functions spherical (lon=x, lat=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  #save declination and inc and calculate new system for rotation
  newSlon <- ifelse((lon+180)>360,lon-180,lon+180)
  newSlat <- 90-lat
  newSlonr <- d2r(newSlon)
  newSlatr <- d2r(newSlat)
  a95 <- A
  circle <- as.data.frame(matrix(ncol=2,nrow=0))
  #loop that create a95 and rotate it around new coordinate (dec, inc)
  for (i in seq(0,360,2)){
    circleP <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(i,(90-a95))
    y <- s2cy(i,(90-a95))
    z <- s2cz(90-a95)
    vec <- as.matrix(c(x,y,z))
    R_elements <- c(cos(newSlatr)*cos(newSlonr), -sin(newSlonr), -sin(newSlatr)*cos(newSlonr),
                    cos(newSlatr)*sin(newSlonr), cos(newSlonr), -sin(newSlatr)*sin(newSlonr),
                    sin(newSlatr), 0, cos(newSlatr))
    R <- matrix(R_elements,nrow=3, byrow=TRUE)
    newvec <- R%*%vec
    newlon <- r2d(atan2(newvec[2,1],newvec[1,1]))
    newlon <- ifelse(newlon<0,newlon+360,newlon)
    #absolute value avoid point outside the graph
    newlat <- r2d(asin(newvec[3,1]))
    circleP[1,1:2] <- c(newlon,newlat)
    circle <- rbind(circle,circleP)
  }
  colnames(circle) <- c("lon","lat")
  circle$x <- c2x(circle$lon,circle$lat)
  circle$y <- c2y(circle$lon,circle$lat)
  circle$cut <- cut(circle$lon,circle$lat)
  #restore screen
  par(fig=c(0,1,0,1))
  #standalone graph or on existing graph
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    sph_ortho(lat = lat0,long = lon0,grid = grid, coast=coast)
    }

  X <- c2x(lon,lat)
  Y <- c2y(lon,lat)
  CUT <- cut(lon,lat)

  #plot alfa 95
  polygon(circle$x,circle$y, col=col_A, lwd=0.8,lty= ifelse(CUT>0,1,3))

  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help (type ?plot_PA95) for info.",call. = F)}

  if(CUT>0){
    points(X,Y, pch=pch,cex=size, col="black",
           bg= col_f)
  }else{
    points(X,Y, pch=pch,cex=size, col="black",
           bg=col_b)
  }
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot great circle on spherical projection with given pole and camera location
plot_plane_sph <- function(P_long=0,P_lat=0,lon0=0,lat0=90,plot_pole=TRUE,on_plot=TRUE,col_f="red",col_b="white",coast=F){
  #service functions
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #plot empty sph if required
  if(on_plot==FALSE){PmagDiR::sph_ortho(lat = lat0,long = lon0,coast = coast)}

  #create equatorial circle
  eq_circle <- data.frame(matrix(ncol=2,nrow = 181))
  colnames(eq_circle) <- c("lon","lat")
  eq_circle[1:181,1] <- seq(0,360,2)
  eq_circle[1:181,2] <- 0
  #rotate circle to real coordinates
  r_circle <- PmagDiR::bed_DI(DI = eq_circle,in_file = F,bed_az = (P_long+180)%%360,bed_plunge = (P_lat-90),export = F)
  colnames(r_circle) <- c("lon","lat")
  #transform coordinates of circle in x and y file
  r_circle$x <- c2x(lon = r_circle[,1],lat = r_circle[,2])
  r_circle$y <- c2y(lon = r_circle[,1],lat = r_circle[,2])
  r_circle$cut <- cut(lon = r_circle[,1],lat = r_circle[,2])
  r_circle <- r_circle[,-c(1,2)]
  l <- 1
  i <- 1
  #double loop to that breaks table depending on sign cut, to avoid line accross the globe
  repeat{
    repeat{
      if(sign(r_circle[l+1,3])!=sign(r_circle[l,3]) || l==nrow(r_circle)){
        if(sign(r_circle[l,3])>=0) {points(x = r_circle[i:l,1],y = r_circle[i:l,2],type="l", col="blue")}
        else if(sign(r_circle[l,3])<0) {points(x = r_circle[i:l,1],y = r_circle[i:l,2],type="l", lty=2,col="blue")}
        break
      }
      l <- l+1
    }
    if(l==nrow(r_circle)) break
    i <- l
    l <- l+1
  }
  if(plot_pole==T){PmagDiR::plot_PA95(lon = P_long,lat = P_lat,A = 0,
                                      lon0 = lon0,lat0 = lat0,col_f = col_f,col_b = col_b,
                                      symbol = "s",on_plot = TRUE)}
}

#plot pole with A95 and Apparent polar wander path
plot_pole_APWP <- function(lon,lat,A,lon0=0,lat0=90,grid=30, col_f="red",col_b="white",col_l="black",col_A=rgb(1,0,0,0.30), symbol="c",coast=FALSE, on_plot=FALSE, save=FALSE, name="A95",APWP="V23", S_APWP=FALSE){
  #plot pole
  plot_PA95(lon=lon, lat = lat,A = A,lon0=lon0, lat0=lat0, grid=grid, col_f = col_f, col_b= col_b, col_l=col_l, col_A = col_A, symbol=symbol, coast=coast, on_plot = on_plot,save=FALSE)

  #plot APWP if requested during process
  pAPWP <- readline("Plot APWP? (y or n): ")
  if(pAPWP=="y"){
    plot_APWP(APWP=APWP, lon0 = lon0,lat0 = lat0,grid = grid,on_plot = T,S_APWP = S_APWP)
  }
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot only apparent polar wander path
plot_APWP <- function(APWP= "V23",lon0=0,lat0=90,grid=30,col="gray",symbol="c",size=0.6, coast=FALSE, on_plot=FALSE, save=FALSE, name="APWP",S_APWP=FALSE,Shiny=FALSE,Y=0,O=320,frame=1,Age_size=1){
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    sph_ortho(lat = lat0,long = lon0,grid = grid, coast=coast)
  }

  #plot APWP
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #questions on APW age and frame
  if(Shiny==FALSE){
    if(APWP=="V23"){
      cat("Frames:
(1) South Africa
(2) North America
(3) South America
(4) Europe
(5) India
(6) Australia
(7) Antarctica
(8) Pacific (0 to 80_Ma)
(9) Iberia (0 to 80 Ma)")
    } else if (APWP=="T12"){
      cat("Frames:
(1) South Africa
(2) North America
(3) Europe
(4) India
(5) Amazonia
(6) Australia
(7) East Antarctica")
    }
    frame <- as.numeric(readline("insert frame (number): "))
    cat("APWP range from 0 to 320 Ma every 10 Myr.
")
    Y <- round(as.numeric(readline("Insert younger age: ")),-1)
    O <- round(as.numeric(readline("Older age: ")),-1)
  }
  if(Shiny==TRUE){
    Y=Y
    O=O
    frame=frame
  }
  col1 <- (frame*2)+1
  col2 <- (frame*2)+2
  if(is.na(O)==TRUE){O <- 320}
  if(O>320){O <- 320}
  if(is.na(Y)==TRUE){Y <- 0}
  #if frame is 8 or 9 (only in V23) and age is too old it fixes it
  if(frame==8 && O>80) O <- 80
  if(frame==9 && O>80) O <- 80
  Y <- (Y/10)+1
  O <- (O/10)+1
  #select apwp file
  if(APWP=="V23") G <- V23_GAPWP
  if(APWP=="T12") G <- T12_GAPWP
  #flip if necessary
  if(S_APWP==FALSE) {G[,col1:col2] <- flip_DI(G[,col1:col2])}
  par(fig=c(0,1,0,1), new=TRUE)
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  #line connecting APWP
  lin <- as.data.frame(c2x(G[Y:O,col1],G[Y:O,col2]))
  colnames(lin) <- "lx"
  lin$ly <- c2y(G[Y:O,col1],G[Y:O,col2])
  lin$cut <- cut(G[Y:O,col1],G[Y:O,col2])
  lines(lin$lx,lin$ly,cex=1)
  #plot poles APWP
  for (i in Y:O){
    plot_PA95(lon = G[i,col1],lat = G[i,col2],A = G[i,2],lon0 = lon0,lat0 = lat0,on_plot = T,col_f = col,symbol=symbol,size=size, col_l = "black",col_A=rgb(1,0.9,0,0.30))
  }
  text1 <- paste(G[Y,1],"Ma")
  text2 <- paste(G[O,1], "Ma")
  text(x=lin[1,1], y=lin[1,2],pos=4,substitute(paste(bold(text1))), cex= Age_size)
  text(x=lin[length(lin$lx),1], y=lin[length(lin$lx),2],pos=4,substitute(paste(bold(text2))), cex= Age_size)

  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function plotting directions
plot_DI <- function(DI,single_mode=FALSE, down=TRUE,symbol="c", col_d="blue",col_u="cyan",col_ext="black", on_plot=FALSE, title="",save=FALSE,name="Equal_area"){
  library("dplyr", warn.conflicts = FALSE)

  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}

  data <- DI[,1:2]
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec","inc")
  if(single_mode==TRUE) {data <- common_DI(data,down=down)}
  data_U <- filter_all(data,all_vars(inc<0))
  data_D <- filter_all(data,all_vars(inc>=0))
  data_U$inc <- abs(data_U$inc)
  xU <- a2cx(data_U$inc,data_U$dec)
  yU <- a2cy(data_U$inc,data_U$dec)
  xD <- a2cx(data_D$inc,data_D$dec)
  yD <- a2cy(data_D$inc,data_D$dec)
  if(on_plot==FALSE){
    equalarea(title=title)
  }
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  points(xD,yD, pch=pch,col=col_ext,
         bg= col_d)
  points(xU,yU, pch=pch,col=col_ext,
         bg=col_u)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot plain intersection in equal area given declination and inclination of pole
plot_plane <- function(D,I, col_cD="black",col_cU="grey", pole=TRUE, col_d="red",col_u="white", symbol="s",on_plot=TRUE,save=FALSE,name="plane"){
  #functions degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #functions spherical (Dec=x, Inc=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  #set intial declination of circle, to avoid lines when plotted
  #strike of pole dec
  data <- D-90
  data <- as.data.frame(seq(data,data+359,1))
  data$inc <- 0
  colnames(data) <- c("dec","inc")
  #define plane Azimuth and dip from pole
  Paz <- D+180
  Pdip <- 90-I
  #sines and cosines of plane coord
  sbd <- -sin(d2r(Paz))
  cbd <- cos(d2r(Paz))
  sbi <- sin(d2r(Pdip))
  cbi <- cos(d2r(Pdip))
  #create new rotated circle
  newDI <- as.data.frame(matrix(ncol=2,nrow=0))
  for(i in 1:length(data[,1])){
    newDI_p <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(data[i,1],data[i,2])
    y <- s2cy(data[i,1],data[i,2])
    z <- s2cz(data[i,2])
    xn <- x*(sbd^2+cbd^2*cbi)+
      y*(cbd*sbd*(1-cbi))+
      z*sbi*cbd
    yn <- x*cbd*sbd*(1-cbi)+
      y*(cbd^2+sbd*sbd*cbi)-
      z*sbd*sbi
    zn <- -(x*cbd*sbi-
              y*sbi*sbd-
              z*cbi)
    newdec <- r2d(atan2(yn,xn))
    newdec <- ifelse(newdec<0,newdec+360,newdec)
    newinc <- r2d(asin(zn))
    newDI_p[1,1:2] <- c(newdec,newinc)
    newDI <- rbind(newDI,newDI_p)
  }
  #convert inc to absolute
  newDI[,2] <- abs(newDI[,2])
  #positive inclination half circle
  circle_U <- newDI[1:180,]
  #negative inclination half circle
  circle_D <- newDI[181:360,]
  colnames(circle_U) <- c("dec","inc")
  colnames(circle_D) <- c("dec","inc")
  #convert to x y
  circle_U$x <- a2cx(circle_U$inc,circle_U$dec)
  circle_U$y <- a2cy(circle_U$inc,circle_U$dec)
  circle_D$x <- a2cx(circle_D$inc,circle_D$dec)
  circle_D$y <- a2cy(circle_D$inc,circle_D$dec)
  #standalone graph or on existing graph
  if (on_plot==FALSE) equalarea()
  UD <- ifelse(I>0,"D","U")
  I <- abs(I)
  X <- a2cx(I,D)
  Y <- a2cy(I,D)
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Type ?plot_DI for info.",call. = F)}

  #plot pole only if pole==TRUE
  if(pole==TRUE){
    if(UD=="D"){
      points(X,Y, pch=pch,cex=1.3, col="black",
             bg= col_d)
    }else{
      points(X,Y, pch=pch,cex=1.3, col="black",
             bg=col_u)
    }
  }
  points(x=circle_U$x,y=circle_U$y,type="l", col=col_cU,lty=2)
  points(x=circle_D$x,y=circle_D$y,type="l", col=col_cD)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot virtual geomagnetic poles
plot_VGP <- function(VGP,lat=90,long=0,grid=30, col="black", on_plot=FALSE,auto_cent=TRUE,exp=TRUE,coast=FALSE, title="",save=TRUE,A95=FALSE,name="VGP"){
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #manipulate data
  #VGP <- na.omit(VGP)
  colnames(VGP) <- c("lon","lat")
  vgpsN <- common_DI(VGP,down = ifelse(mean(VGP$lat)<0,FALSE,TRUE))
  PPole <- fisher(vgpsN)

  #fix point of view
  if(auto_cent==FALSE){
    #center of proj is Lon0 & Lat0
    lon0 <- long
    lat0 <- lat
  }else{
    lon0 <- PPole[1,1]
    lat0 <- PPole[1,2]
  }
  if(on_plot==FALSE){
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    sph_ortho(lat=lat0,long=lon0,grid=grid,coast=coast, title=title)
    }

  coord <- as.data.frame(lon0)
  coord$lat0 <- lat0
  if(exp==TRUE){
    write.csv(round(coord, digits=2),file="center_coordinates.csv", row.names = F)
  }
  cat(paste("Center coordinates:

"))
  print(round(coord,digits=2),row.names=F)

  VGP$x <- c2x(VGP$lon,VGP$lat)
  VGP$y <- c2y(VGP$lon,VGP$lat)
  VGP$cut <- cut(VGP$lon,VGP$lat)

  points(VGP$x,VGP$y,pch=ifelse(VGP$cut>0,21,1),col="black",bg=col)
  if(A95==TRUE){
    plot_PA95(lon = PPole[1,1],lat = PPole[1,2],A = PPole[1,3],lon0 = lon0,lat0 = lat0,on_plot = TRUE,symbol = "d",col_l = "red")
    text <- paste("N: ",PPole[1,4],"
Long: ", round(PPole[1,1],digits=2),"
Lat: ", round(PPole[1,2], digits=2),"
A95: ", round(PPole[1,3], digits=2))
    text(x=0.75, y=-0.85,pos=4,text, cex= 0.85)
  }
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#reversal test boostrapped following Tauxe
revtest <- function(DI,nb=1000,export=TRUE, name="reversal_test"){
  #fucnctions deg to rads and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")

  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))

  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)

  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))

  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)

  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T,symmetric = TRUE)
  T_vec <- T_e$vectors
  T_val <- T_e$values

  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)

  #flip V1 if negative
  V1dec <- ifelse(V1inc<0,ifelse((V1dec+180)>360,V1dec-180,V1dec+180),V1dec)
  V1inc <- ifelse(V1inc<0,-V1inc,V1inc)


  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  m1ind <- as.numeric(which(data$diff<=90), arr.ind = TRUE)
  m2ind <- as.numeric(which(data$diff>90), arr.ind = TRUE)

  #terminate if distribution is not bimodal
  if(length(m2ind)<1) stop("
DISTRIBUTION NOT BIMODAL")
  mode1 <- data[m1ind,1:2]
  mode2 <- data[m2ind,1:2]

  #flip mode 2 same as mode 1
  mode2$dec <- ifelse((mode2$dec+180)>360,mode2$dec-180,mode2$dec+180)
  mode2$inc <- -mode2$inc
  nb <- nb
  n <- 0
  mode1B <- as.data.frame(matrix(ncol=3, nrow=0))
  mode2B <- as.data.frame(matrix(ncol=3, nrow=0))

  #simulate pseudosamples of mode 1
  repeat{
    n <- n+1
    mode1B_p <- as.data.frame(matrix(ncol=3, nrow=1))
    Bdata <- boots_DI(mode1)
    Bdata$x <- cos(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$y <- sin(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$z <- sin(d2r(Bdata$inc))
    mode1B_p[1,1] <- mean(Bdata$x)
    mode1B_p[1,2] <- mean(Bdata$y)
    mode1B_p[1,3] <- mean(Bdata$z)
    mode1B <- rbind(mode1B,mode1B_p)
    if(((n%%50)==0)==TRUE){
      cat(paste(n,"simulations out of",nb,"of mode 1 done
"))
    }
    if(n==nb) break
  }
  colnames(mode1B) <- c("x","y","z")
  mode1B$dec <- r2d(atan2(mode1B$y,mode1B$x))
  mode1B$dec <- mode1B$dec%%360
  mode1B$inc <- r2d(asin(mode1B$z))
  n <- 0
  #simulate pseudosamples of mode 2
  repeat{
    n <- n+1
    mode2B_p <- as.data.frame(matrix(ncol=3, nrow=1))
    Bdata <- boots_DI(mode2)
    Bdata$x <- cos(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$y <- sin(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$z <- sin(d2r(Bdata$inc))
    mode2B_p[1,1] <- mean(Bdata$x)
    mode2B_p[1,2] <- mean(Bdata$y)
    mode2B_p[1,3] <- mean(Bdata$z)
    mode2B <- rbind(mode2B,mode2B_p)
    if(((n%%50)==0)==TRUE){
      cat(paste(n,"simulations out of",nb,"of mode 2 done
"))
    }
    if(n==nb) break
  }
  colnames(mode2B) <- c("x","y","z")
  mode2B$dec <- r2d(atan2(mode2B$y,mode2B$x))
  mode2B$dec <- ifelse(mode2B$dec<0,mode2B$dec+360,mode2B$dec)
  mode2B$inc <- r2d(asin(mode2B$z))

  #isolate components of models
  B1x <- sort(mode1B[,1])
  B1y <- sort(mode1B[,2])
  B1z <- sort(mode1B[,3])
  B2x <- sort(mode2B[,1])
  B2y <- sort(mode2B[,2])
  B2z <- sort(mode2B[,3])

  #define low and high boostrapped margins
  confn <- 0.95
  num <- round((nb*(1-confn))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  B1x_l <- c(B1x[Lconf],B1x[Uconf])
  B2x_l <- c(B2x[Lconf],B2x[Uconf])
  B1y_l <- c(B1y[Lconf],B1y[Uconf])
  B2y_l <- c(B2y[Lconf],B2y[Uconf])
  B1z_l <- c(B1z[Lconf],B1z[Uconf])
  B2z_l <- c(B2z[Lconf],B2z[Uconf])

  #max and min values for graphs
  xmax <- round(max(c(B1x,B2x)), digits=1)+0.05
  xmin <- round(min(c(B1x,B2x)),digits=1)-0.05
  ymax <- round(max(c(B1y,B2y)), digits=1)+0.05
  ymin <- round(min(c(B1y,B2y)),digits=1)-0.05
  zmax <- round(max(c(B1z,B2z)), digits=1)+0.05
  zmin <- round(min(c(B1z,B2z)),digits=1)-0.05

  #function that extract intervals and counts from hist function and make cumulative curve
  cumulative_curve <- function(x){
    h <- hist(x, breaks=50,plot = FALSE)
    cnts <- h[["counts"]]
    t <- length(cnts)
    new_c <- as.data.frame(matrix(ncol=1,nrow = 1))
    for(i in 1:t){
      if(i==1){new_cp <- cnts[1]}
      if(i>1) {new_cp <- new_c[i,1]+cnts[i-1]}
      new_c <- rbind(new_c,new_cp)
    }
    new_c <- na.omit(new_c)
    breaks <- as.data.frame(h[["mids"]])
    cumul <- cbind(breaks,new_c)
    colnames(cumul) <- c("breaks","counts")
    cumul$counts <- cumul$counts/nb
    return(cumul)
  }
  cu1x <- cumulative_curve(B1x)
  cu2x <- cumulative_curve(B2x)
  cu1y <- cumulative_curve(B1y)
  cu2y <- cumulative_curve(B2y)
  cu1z <- cumulative_curve(B1z)
  cu2z <- cumulative_curve(B2z)
  text1 <- "Equal area projections"
  text2 <- "Normalized cumulative distributions"

  #clean screen to avoid figure over figure
  par(fig=c(0,1,0,1))
  plot(0, xlim=c(0,1), ylim=c(0,1),
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

  #plot title for equal area
  par(fig=c(0,0.65,0.4,1))
  plot(NA, xlim=c(0,1), ylim=c(0,1),
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  text(x=0.5,y=1,text1,cex=1.2)

  #plot title for cumulative distributions
  par(fig=c(0.55,1,0.5,1),new=TRUE)
  plot(NA, xlim=c(0,1), ylim=c(0,1),
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  text(x=0.5,y=0.9,text2,cex=1.2)

  #plot equal area
  #real data
  par(fig=c(0,0.65,0.4,0.98),new=TRUE)
  plot_DI(mode1[,1:2],col_d = rgb(1,0,0,0.30),col_u=rgb(1,0.75,1,0.30),col_ext = NA, title = "Data common mode")
  plot_DI(mode2[,1:2],col_d = rgb(0,0,1,0.30),col_u=rgb(0,1,1,0.30),col_ext = NA,on_plot = TRUE)
  #pseudosamples mean
  par(fig=c(0,0.65,0,0.58), new=TRUE)
  plot_DI(mode1B[,4:5],col_d = rgb(1,0,0,0.20),col_u=rgb(1,0.75,1,0.30),col_ext = NA,title = "Pseudosample means")
  plot_DI(mode2B[,4:5],col_d =rgb(0,0,1,0.20),col_u=rgb(0,1,1,0.30),col_ext = NA,on_plot = TRUE)

  #plot cumulative distributions
  par(fig=c(0.55,1,0.58,0.91),new=TRUE)
  #plot x
  plot(0,type="n",xlim=c(xmin,xmax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
  title(xlab = "x axis", line=2)
  rect(xleft = B1x_l[1],xright=B1x_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
  rect(xleft = B2x_l[1],xright=B2x_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
  points(cu1x, type="l", lwd=2, col="red")
  points(cu2x, type="l", lwd=2, col="blue")

  par(fig=c(0.55,1,0.34,0.67),new=TRUE)
  #plot y
  plot(0,type="n",xlim=c(ymin,ymax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
  title(xlab = "y axis", line=2)
  rect(xleft = B1y_l[1],xright=B1y_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
  rect(xleft = B2y_l[1],xright=B2y_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
  points(cu1y, type="l", lwd=2, col="red")
  points(cu2y, type="l", lwd=2, col="blue")

  par(fig=c(0.55,1,0.10,0.43),new=TRUE)
  #plot z
  plot(0,type="n",xlim=c(zmin,zmax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
  title(xlab = "z axis", line=2)
  rect(xleft = B1z_l[1],xright=B1z_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
  rect(xleft = B2z_l[1],xright=B2z_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
  points(cu1z, type="l", lwd=2, col="red")
  points(cu2z, type="l", lwd=2, col="blue")
  #reset screen
  par(fig=c(0,1,0,1))
  #export if requested
  if(export==TRUE){
    save_pdf(name=paste(name,".pdf"),height =8,width =11 )
    boot_stat <- matrix(c(B1x_l[1],B1x_l[2],B2x_l[1],B2x_l[2],
                          B1y_l[1],B1y_l[2],B2y_l[1],B2y_l[2],
                          B1z_l[1],B1z_l[2],B2z_l[1],B2z_l[2]),nrow = 3,byrow = TRUE)
    rownames(boot_stat) <- c("x","y","z")
    colnames(boot_stat) <- c("mode1_L","mode1_H","mode2_L","mode2_H")
    write.csv(round(boot_stat,digits = 3),paste(name,"bootstrap_stat.csv"),row.names = TRUE)
    cat("Figure saved as ",paste(name,".pdf"),"
")
  }
}

#function that rotate Longitude and Latitude around a Euler pole
rot_DI <- function(Lonlat,P_long=0,P_lat=90,rot=0){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions spherical (long=x, lat=y) to Cartesian
  s2c1 <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2c2 <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2c3 <- function(y) {sin(d2r(y))}

  #functions converting cartesian to spherical
  c2sLon <- function(C1,C2) {r2d(atan2(C2,C1))}
  c2sLat <- function(C1,C2,C3) {r2d(asin(C3/(sqrt((C1^2)+(C2^2)+(C3^2)))))}

  #define cos and sin of rotation for simplicity
  cosrot <- cos(d2r(rot))
  sinrot <- sin(d2r(rot))
  #build rotation matrix
  R <- matrix(c((s2c1(P_long,P_lat)^2)*(1-cosrot)+cosrot,
                s2c1(P_long,P_lat)*s2c2(P_long,P_lat)*(1-cosrot)-s2c3(P_lat)*sinrot,
                s2c1(P_long,P_lat)*s2c3(P_lat)*(1-cosrot)+s2c2(P_long,P_lat)*sinrot,
                s2c2(P_long,P_lat)*s2c1(P_long,P_lat)*(1-cosrot)+s2c3(P_lat)*sinrot,
                (s2c2(P_long,P_lat)^2)*(1-cosrot)+cosrot,
                s2c2(P_long,P_lat)*s2c3(P_lat)*(1-cosrot)-s2c1(P_long,P_lat)*sinrot,
                s2c3(P_lat)*s2c1(P_long,P_lat)*(1-cosrot)-s2c2(P_long,P_lat)*sinrot,
                s2c3(P_lat)*s2c2(P_long,P_lat)*(1-cosrot)+s2c1(P_long,P_lat)*sinrot,
                (s2c3(P_lat)^2)*(1-cosrot)+cosrot),
              nrow = 3,ncol = 3,byrow = F)
  #creates result file
  Lonlat_R <- data.frame(matrix(ncol = 2,nrow = 0))
  #apply rotation to all data
  for(i in 1:nrow(Lonlat)){
    C <- matrix(c(s2c1(Lonlat[i,1],Lonlat[i,2]),
                  s2c2(Lonlat[i,1],Lonlat[i,2]),
                  s2c3(Lonlat[i,2])),
                nrow = 3,ncol = 1)
    C_rot <- R%*%C
    Lonlat_R_temp <- data.frame(t(c(c2sLon(C_rot[1,1],C_rot[2,1])%%360,
                                    c2sLat(C_rot[1,1],C_rot[2,1],C_rot[3,1]))))
    Lonlat_R <- rbind(Lonlat_R,Lonlat_R_temp)
  }
  colnames(Lonlat_R) <- c("Long_R","Lat_R")
  return(Lonlat_R)
}

#pdf printing standard size
save_pdf <- function(name="Figure.pdf",width=11,height=8){
  dev.print(pdf,name,width = width, height = height)
}

#plot spherical ortographic projection centered in specified coordinates
sph_ortho <- function(lat=90,long=0,grid=30,coast=FALSE, title="") {
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  #center of proj is Lon0 & Lat0
  lon0 <- long
  lat0 <- lat
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}

  #fix frame
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  #plot coastline if true
  if(coast==TRUE){
    cst <- world_coastline
    colnames(cst) <- c("lon","lat")
    cst$x <- ifelse(cut(cst$lon,cst$lat)<0,NA,c2x(cst$lon,cst$lat))
    cst$y <- ifelse(cut(cst$lon,cst$lat)<0,NA,c2y(cst$lon,cst$lat))
    lines(cst$x,cst$y,col="black",lwd=0.75)
  }
  #if grid==0 does not plot parallels and meridians
  if(grid!=0){
    #plot_main_parallel
    #longitude circle
    lats <- seq(-(90-grid),(90-grid),grid)
    for(i in lats){
      lon_lat_p <-  as.data.frame(0:360)
      lon_lat_p$lat <- rep(i)
      lon_lat_p$x <- ifelse(cut(lon_lat_p[,1],lon_lat_p[,2])<0,NA,
                            c2x(lon_lat_p[,1],lon_lat_p[,2]))
      lon_lat_p$y <- ifelse(cut(lon_lat_p[,1],lon_lat_p[,2])<0,NA,
                            c2y(lon_lat_p[,1],lon_lat_p[,2]))
      lines(lon_lat_p$x,lon_lat_p$y,col="gray", pch=16, cex=0.3)
    }

    #plot_main_meridians
    lons <- seq(grid,360,grid)
    for(i in lons){
      lat_lon_m <- as.data.frame(seq(-89,89,1))
      lat_lon_m$lon <- rep(i)
      lat_lon_m$x <- ifelse(cut(lat_lon_m[,2],lat_lon_m[,1])<0,NA,
                            c2x(lat_lon_m[,2],lat_lon_m[,1]))
      lat_lon_m$y <- ifelse(cut(lat_lon_m[,2],lat_lon_m[,1])<0,NA,
                            c2y(lat_lon_m[,2],lat_lon_m[,1]))
      lines(lat_lon_m$x,lat_lon_m$y,col="gray", pch=16, cex=0.3)
    }
  }

  #plot black frame around globe
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #plot external circle
  frame_dec = 0:360
  frame_inc=rep(0,length(frame_dec))
  x = a2cx(frame_inc,frame_dec)
  y = a2cy(frame_inc,frame_dec)
  lines(x, y, col = "black")
  title(xlab = title, line=0.2, cex=0.1)
}

#function that deform dec Inc from eigenvectors
strain_DI <- function(DIAP,M,export=FALSE,name="strained_dirs"){
  library(matlib)
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #load directions data
  data <- DIAP
  data <- na.omit(data)
  dirs <- data[,1:2]
  bed <- data[,3:4]
  colnames(dirs) <- c("dec", "inc")
  colnames(bed) <- c("b_az","b_plunge")

  #directions in Cartesian coordinates
  dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$z <- sin(d2r(dirs$inc))

  #bedding in cartesian coordinates
  bed$x <- cos(d2r(bed$b_az))*cos(d2r(bed$b_plunge))
  bed$y <- sin(d2r(bed$b_az))*cos(d2r(bed$b_plunge))
  bed$z <- sin(d2r(bed$b_plunge))

  new_DI <- data.frame(matrix(ncol=2,nrow=0))
  new_bed <- data.frame(matrix(ncol=2,nrow=0))

  for (i in 1:nrow(data)){
    #strain dirs
    dircart <- t(as.matrix(dirs[i,3:5]))
    strain <- M%*%dircart
    NewInc <- r2d(asin(strain[3,1]/(sqrt((strain[1,1]^2)+(strain[2,1]^2)+(strain[3,1]^2)))))
    NewDec <- r2d(atan2(strain[2,1],strain[1,1]))
    NewDec <- NewDec%%360

    new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc))
    new_DI <- rbind(new_DI,new_DI_p)

    #strain bed
    bedcart <- t(as.matrix(bed[i,3:5]))
    b_strain <- M%*%bedcart
    New_pl <- r2d(asin(b_strain[3,1]/(sqrt((b_strain[1,1]^2)+(b_strain[2,1]^2)+(b_strain[3,1]^2)))))
    New_az <- r2d(atan2(b_strain[2,1],b_strain[1,1]))
    New_az <- New_az%%360

    new_bed_p <- cbind(as.data.frame(New_az),as.data.frame(New_pl))
    new_bed <- rbind(new_bed,new_bed_p)
  }
  colnames(new_DI) <- c("str_dec","str_inc")
  colnames(new_bed) <- c("str_B_az","str_Binc")
  str_data <- cbind(new_DI,new_bed)
  if(export==TRUE){write.csv(round(str_data,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(str_data)
}

#function generating E from TK03.GAD model from given inclination
tk03 <- function(I) {
  x <- sqrt(I^2)
  y <- 2.895-(0.01466*x)-(0.0003525*(x^2))+(0.00000316*(x^3)) #E-I equation
  return(y)
}

#unflat directions with given f
unflat_DI <- function(DI,f,export=FALSE,name="unflattened_dirs.csv") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  DI[,2]<- r2d(atan(tan(d2r(DI[,2]))/f))
  if(export==TRUE){write.csv(round(DI,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(DI)
}

#unstrain bootstrpped pseudosamples of directions
unstr_boot <- function(unstr_file,nb= 100,S_vec,Lin,Fol,ns=1,confidence=95,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE,hist=TRUE,save=TRUE,name="Unstrain_bootstrap"){
  cat("
!!Unstrain of pseudosamples is SLOW!!

")
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #load directions data
  dat <- unstr_file[[1]]
  dat <- na.omit(dat)
  colnames(dat) <- c("dec", "inc","baz","binc")

  #load inc_e_dec result
  inc_e_dec <- unstr_file[[5]]

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(inc_e_dec[1,2]>3.5, ceiling(inc_e_dec[1,2]),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="red", lwd=3)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  points(x=inc_e_dec[length(inc_e_dec[,1]),1],y=inc_e_dec[length(inc_e_dec[,1]),2],
         pch=21, col="black", bg="red", cex=1.2)

  #text for figure
  N <- length(dat[,1])
  Inc <- round(inc_e_dec[1,1],digits = 1)
  Ecut <- round(inc_e_dec[1,2],digits = 2)
  Edec <- round(inc_e_dec[1,3],digits = 1)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",Edec)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #unstrained parameters for text2 of figure
  inc_nstr <- round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1)
  E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
  Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

  #put text unstrained in figure
  text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  #create empty files
  init_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  final_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  colnames(init_E_I) <- c("Inc","E","E_dec")
  colnames(final_E_I)<- c("Inc","E","E_dec")

  #start bootstrapping
  n <- 0
  repeat{
    n <- n+1
    Seq_I_E_b <- as.data.frame(matrix(ncol=3,nrow=0))
    E_declin <- as.data.frame(matrix(ncol=1,nrow=0))
    data <- boots_DI(dat)
    dirs <- data[,1:2]
    bed <- data[,3:4]

    #directions in Cartesian coordinates
    dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
    dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
    dirs$z <- sin(d2r(dirs$inc))

    #bedding in Cartesian coordinates
    bed$x <- cos(d2r(bed$baz))*cos(d2r(bed$binc))
    bed$y <- sin(d2r(bed$baz))*cos(d2r(bed$binc))
    bed$z <- sin(d2r(bed$binc))

    #set parameters of deforming matrix
    Lincr <- (Lin-1)/ns
    Fincr <- (Fol-1)/ns
    L <- 1
    F <- 1

    #gradual unstrain of the n pseudosample
    repeat{
      #anisotropy degree
      P <- F*L
      #eigenvalues
      K1 <- (3*L)/(L+(1/F)+1)
      K2 <- K1/L
      K3 <- K1/P

      #matrix of new eigenvalue
      M <- c(K1,0,0,0,K2,0,0,0,K3)
      M <- matrix(M,nrow=3,byrow=T)

      #combines given eigenvalues with Strain directions
      S <- S_vec%*%M%*%inv(S_vec)
      new_DI <- data.frame(matrix(ncol=4,nrow=0))
      for (i in 1:length(data[,1])){
        #unstrain dirs
        dircart <- t(as.matrix(dirs[i,3:5]))
        unstr <- S%*%dircart
        NewInc <- r2d(asin(unstr[3,1]/(sqrt((unstr[1,1]^2)+(unstr[2,1]^2)+(unstr[3,1]^2)))))
        NewDec <- r2d(atan2(unstr[2,1],unstr[1,1]))
        NewDec <- ifelse(NewDec<0,NewDec+360,NewDec)

        #unstrain bedding
        bedcart <- t(as.matrix(bed[1,3:5]))
        bunstr <- S%*%bedcart
        Newbinc <- r2d(asin(bunstr[3,1]/(sqrt((bunstr[1,1]^2)+(bunstr[2,1]^2)+(bunstr[3,1]^2)))))
        Newbaz <- r2d(atan2(bunstr[2,1],bunstr[1,1]))
        Newbaz <- ifelse(Newbaz<0,Newbaz+360,Newbaz)

        new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc),
                          as.data.frame(Newbaz),as.data.frame(Newbinc))
        new_DI <- rbind(new_DI,new_DI_p)
      }
      colnames(new_DI) <- c("dec","inc","baz","binc")
      NBdecInc <- bed_DI(new_DI)
      inc_e_dec_p <- inc_E_finder(NBdecInc)
      Seq_I_E_b <- rbind(Seq_I_E_b,inc_e_dec_p)
      #take only absolute value of declination for the breaking conditions
      E_declin <- rbind(E_declin,abs(inc_e_dec_p[1,3]))

      #if cross== TRUE breaks loop when tk03.GAD is crossed
      if(cross==TRUE)  {if(any(Seq_I_E_b$E>tk03(Seq_I_E_b$V1inc))
                           && any(Seq_I_E_b$E<tk03(Seq_I_E_b$V1inc))) break}

      #if EdMAX==TRUE breaks loop when max of Edec is reached
      if(EdMAX==TRUE && length(E_declin[,1])>3){
        l <- length(E_declin[,1])
        if(E_declin[l,1]<E_declin[l-1,1] && E_declin[l-1,1]>E_declin[l-2,1]){    ##HERE
          E_declin <- E_declin[-l,]
          break
        }
      }
      #break loops when minimu Edec reached
      if(EdMIN==TRUE && length(E_declin[,1])>3){
        l <- length(E_declin[,1])
        if(E_declin[l,1]>E_declin[l-1,1] && E_declin[l-1,1]<E_declin[l-2,1]){
          E_declin <- E_declin[-l,]
          break
        }
      }

      #break loop if reaches defined F and L
      if(round(L,digits = 4) == round(Lin,digits = 4) &&
         round(F,digits = 4) == round(Fol,digits = 4)) break
      #lineation
      L <- L+Lincr
      #foliation
      F <- F+Fincr
    }

    #if cross==TRUE select only curves that cross tk03.GAD
    if(cross==TRUE){
      if(any(Seq_I_E_b$E>tk03(Seq_I_E_b$V1inc)) && any(Seq_I_E_b$E<tk03(Seq_I_E_b$V1inc))==TRUE){
        points(x=Seq_I_E_b$V1inc, y= Seq_I_E_b$E,
               type= "l", col=rgb(1, 0, 0, 0.15), lwd=1)
        i_E_I <- Seq_I_E_b[1,]
        f_E_I <- Seq_I_E_b[length(Seq_I_E_b[,1]),]
      }
      #if cross== FALSE, use any other condition to fill files
    }else{
      points(x=Seq_I_E_b$V1inc, y= Seq_I_E_b$E,
             type= "l", col=rgb(0, 0.5, 1, 0.15), lwd=1)
      i_E_I <- Seq_I_E_b[1,]
      f_E_I <- Seq_I_E_b[length(Seq_I_E_b[,1]),]
    }
    colnames(i_E_I) <- c("Inc","E","E_dec")
    colnames(f_E_I)<- c("Inc","E","E_dec")
    init_E_I <- rbind(init_E_I,i_E_I)
    final_E_I <- rbind(final_E_I,f_E_I)
    init_E_I <- na.omit(init_E_I)
    final_E_I <- na.omit(final_E_I)

    #message for bootstrapping update
    if(((n%%10)==0)==TRUE) {
      cat(paste(n,"simulations done and",(length(final_E_I[,1])),"pseudosamples saved
"))
    }
    if(length(final_E_I[,1])==nb) {
      cat(paste("Saved",(length(final_E_I[,1])), "pseudosamples after", n,"simulations
"))
      break
    }
  }
  #backup_files
  final_E_Ibk <- final_E_I
  init_E_Ibk <- init_E_I

  #Plot Bootstrapped data
  if(cross==FALSE){
    points(x=final_E_I$Inc,
           y=final_E_I$E,
           pch=16,
           col=rgb(1,0,0,0.7),
           cex=0.5)
  }

  #cut bootstrapped results for 95% (unless different) confidence
  conf <- confidence/100
  num <- round((nb*(1-conf))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  final_E_I <- final_E_I[order(final_E_I$Inc),]
  final_E_I <- final_E_I[Lconf:Uconf,]

  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="yellow", lwd=2)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  #plot cross with EI line
  points(x=inc_e_dec[length(inc_e_dec[,1]),1],y=inc_e_dec[length(inc_e_dec[,1]),2],
         pch=21, col="black", bg=rgb(1,0.4,0.4,1), cex=1.2)

  #draw two lines for 95% confidence margin
  arrows(x0=final_E_I[1,1],x1=final_E_I[1,1],
         y0=1,y1=1.1, lwd=1.5, length = 0)

  arrows(x0=final_E_I[length(final_E_I$Inc),1],
         x1=final_E_I[length(final_E_I$Inc),1],
         y0= 1, y1= 1.1,lwd=1.5, length = 0)

  Inc_l <- round(final_E_I[1,1], digits= 1)
  Inc_u <- round(final_E_I[length(final_E_I$Inc),1], digits=1)

  text(x=final_E_I[1,1], y=1, pos=2, Inc_l)
  text(x=final_E_I[length(final_E_I$Inc),1], y=1, pos=4, Inc_u)

  #plot text results again if covered by red lines
  if(inc_e_dec[1,1]<40){
    text(x=0, y=3.2,pos=4,text, cex= 0.8)
    text(x=20, y=3.2, pos=4, text2, cex=0.8)
  }
  #recalculate boostrtapped confidence for Edec for plotting confidence margin
  final_E_I_Edec <- final_E_Ibk
  colnames(final_E_I_Edec) <- c("Inc","E","E_dec")
  final_E_I_Edec <- final_E_I_Edec[order(final_E_I_Edec$E_dec),]
  final_E_I_Edec <- final_E_I_Edec[Lconf:Uconf,]

  #plot histogram of E_declination with respect V1 before and after correction
  if(hist==TRUE){
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(init_E_Ibk$E_dec, xlim=c(-90,90), breaks= 90,
         axes=FALSE,xlab="",ylab="",col="blue", border="blue", main="")
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(final_E_Ibk$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
         breaks = 90, xlab = "", ylab = "",
         main="", cex.axis=0.8,col="red",border ="red")
    abline(v=final_E_I_Edec[1,3],lwd=1, lty=2)
    abline(v=final_E_I_Edec[length(final_E_I_Edec[,3]),3],lwd=1, lty=2)

    #plot lables closer than standard to axes
    title(xlab = "Edec(°)", line=1.9)
    title(ylab = "Frequency", line=1.9)
  }

  #plot original directions single mode tilt corrected
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(unstr_file[[2]],single_mode = TRUE, title="Original directions")

  #plot unstrained directions tilt corrected
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(unstr_file[[3]], single_mode = TRUE, col_d = "red",col_u = "pink", title = "Unstrained directions")

  #build result file
  stat <- as.data.frame(matrix(ncol=1,nrow=1))
  colnames(stat) <- c("N")
  stat$N <- N
  stat$Inc <- Inc
  stat$E <- Ecut
  stat$Edec <- Edec
  stat$Inc_unstr <- inc_nstr
  stat$Low_inc <- round(final_E_I[1,1], digits = 1)
  stat$Hign_inc <- round(final_E_I[length(final_E_I[,1]),1], digits = 1)
  stat$E_unstr <- E_nstr
  stat$Edec_unstr <- Edec_nstr
  stat$Edec_low <- round(final_E_I_Edec[1,3],digits = 1)
  stat$Edec_high <- round(final_E_I_Edec[length(final_E_I_Edec[,3]),3],digits = 1)

  #save and export results if save=TRUE
  if(save==TRUE){
    cat("
Statistics saved as csv file.
Figure saved as pdf file.

      ")
  }
  #restore screen
  par(fig=c(0,1,0,1))
  if(save==TRUE){
    write.csv(stat,paste(name,"_statistic.csv"), row.names=FALSE)
    save_pdf(paste(name,".pdf"))
  }
}

#file DI require also 3 and 4 columns with Bed_az and bed_plunge, S_vec= strain matrix
unstr_DI <- function(DIAP,S_vec,Lin,Fol,n=1,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE,save=TRUE,name="Unstrain"){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #load directions data
  data <- DIAP
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc","baz","binc")
  dirs <- data[,1:2]
  bed <- data[,3:4]

  #directions in Cartesian coordinates
  dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$z <- sin(d2r(dirs$inc))

  #bedding in Cartesian coordinates
  bed$x <- cos(d2r(bed$baz))*cos(d2r(bed$binc))
  bed$y <- sin(d2r(bed$baz))*cos(d2r(bed$binc))
  bed$z <- sin(d2r(bed$binc))

  inc_e_dec <- as.data.frame(matrix(ncol=3,nrow = 1))
  BdecInc <- bed_DI(data[,1:4])
  inc_e_dec[1,1:3] <- as.data.frame(inc_E_finder(BdecInc))
  inc_e_dec[1,3] <- abs(inc_e_dec[1,3])
  colnames(inc_e_dec) <- c("V1inc","E","DV1V2")

  #set parameters of deforming matrix
  Lincr <- (Lin-1)/n
  Fincr <- (Fol-1)/n
  L <- 1
  F <- 1
  repeat{
    #lineation
    L <- L+Lincr
    #foliation
    F <- F+Fincr
    #anisotropy degree
    P <- F*L
    #eigenvalues
    K1 <- (3*L)/(L+(1/F)+1)
    K2 <- K1/L
    K3 <- K1/P

    #matrix of new eigenvalue
    M <- c(K1,0,0,0,K2,0,0,0,K3)
    M <- matrix(M,nrow=3,byrow=T)

    #combines given eigenvalues with Strain directions
    S <- S_vec%*%M%*%inv(S_vec)
    S_e <- eigen(S,symmetric = T)
    S_vec <- S_e$vectors
    S_val <- S_e$values
    new_DI <- data.frame(matrix(ncol=4,nrow=0))
    for (i in 1:nrow(data)){
      #unstrain dirs
      dircart <- t(as.matrix(dirs[i,3:5]))
      unstr <- S%*%dircart
      NewInc <- r2d(asin(unstr[3,1]/(sqrt((unstr[1,1]^2)+(unstr[2,1]^2)+(unstr[3,1]^2)))))
      NewDec <- r2d(atan2(unstr[2,1],unstr[1,1]))
      NewDec <- NewDec%%360
      #unstrain bedding
      bedcart <- t(as.matrix(bed[1,3:5]))
      bunstr <- S%*%bedcart
      Newbinc <- r2d(asin(bunstr[3,1]/(sqrt((bunstr[1,1]^2)+(bunstr[2,1]^2)+(bunstr[3,1]^2)))))
      Newbaz <- r2d(atan2(bunstr[2,1],bunstr[1,1]))
      Newbaz <- Newbaz%%360

      new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc),
                        as.data.frame(Newbaz),as.data.frame(Newbinc))
      new_DI <- rbind(new_DI,new_DI_p)
    }
    colnames(new_DI) <- c("dec","inc","baz","binc")
    NBdecInc <- bed_DI(new_DI)
    inc_e_dec_p <- inc_E_finder(NBdecInc)
    inc_e_dec <- rbind(inc_e_dec,inc_e_dec_p)
    inc_e_dec[,3] <- abs(inc_e_dec[,3])

    #generate tk03 curve also for later plot
    tkx <- 0:90
    tky <- tk03(tkx)
    if (cross==TRUE){
      #break loop if crosses tk03 line
      if(any(inc_e_dec$E>tk03(inc_e_dec$V1inc)) && any(inc_e_dec$E<tk03(inc_e_dec$V1inc))){
        curve1 <- inc_e_dec[,1:2]
        curve1 <- curve1[order(curve1$E),]
        curve2 <- cbind(as.data.frame(tkx),as.data.frame(tky))
        crxy <- curve_cross(curve1,curve2)
        cat(paste(" Final L: ",round(L,digits = 5),"
", "Final F:",round(F,digits = 5),"
","Final P:",round(F*L,digits = 5)))
        break
      }
    }
    #break loops when maximum Edec reached
    if(EdMAX==TRUE && length(inc_e_dec[,3])>3){
      l <- length(inc_e_dec[,3])
      if(inc_e_dec[l,3]<inc_e_dec[l-1,3] && inc_e_dec[l-1,3]>inc_e_dec[l-2,3]){
        inc_e_dec <- inc_e_dec[-l,]
        break
      }
    }
    #break loops when minimu Edec reached
    if(EdMIN==TRUE && length(inc_e_dec[,3])>3){
      l <- length(inc_e_dec[,3])
      if(inc_e_dec[l,3]>inc_e_dec[l-1,3] && inc_e_dec[l-1,3]<inc_e_dec[l-2,3]){
        inc_e_dec <- inc_e_dec[-l,]
        break
      }
    }
    #break loop if reaches defined F and L
    if(round(L,digits = 4) == round(Lin,digits = 4) &&
       round(F,digits = 4) == round(Fol,digits = 4)) break
  }
  #set figure
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(inc_e_dec[1,2]>3.5, ceiling(inc_e_dec[1,2]),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #text for figure
  N <- as.character(length(data[,1]))
  Inc <- round(inc_e_dec[1,1],digits = 1)
  Ecut <- round(inc_e_dec[1,2],digits = 2)
  V2 <- round(inc_e_dec[1,3],digits = 1)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #plot tk03.GAD model E-I
  points(x=tkx, y= tky, type= "l", col="blue", lwd=3)
  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="red", lwd=3)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  #plot cross with EI line if exists
  if(exists("crxy")==TRUE){
    points(x=crxy[1,1],y=crxy[1,2],
           pch=21, col="black", bg="red", cex=1.2)
  }else{
    points(x=inc_e_dec[length(inc_e_dec[,1]),1],
           y=inc_e_dec[length(inc_e_dec[,1]),2],
           pch=21, col="black", bg="red", cex=1.2)
  }
  #unstrained parameters for text2 of figure
  inc_nstr <- ifelse(exists("crxy")==TRUE,
                     round(crxy[1,1], digits=1),
                     round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1))

  E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
  Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

  #put text unstrained in figure
  text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  #plot Edec during unstrain
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
  plot(x=inc_e_dec$V1inc, y=inc_e_dec$DV1V2, type="l", tck=-0.05,
       frame.plot=FALSE,cex.axis=0.8, mgp=c(1.5,0.5,0.1),
       xlab="",ylab="",yaxt="n",
       ylim=c(floor(min(inc_e_dec$DV1V2)),
                    ceiling(max(inc_e_dec$DV1V2))),
       xlim=c(floor(min(inc_e_dec$V1inc)),
              ceiling(max(inc_e_dec$V1inc))))
  axis(2, cex.axis=0.8, las=2)
  title(xlab = "Inc (°)", line=1.6)
  title(ylab = "Edec (°)", line=2.5)

  points(x=inc_e_dec[1,1],y=inc_e_dec[1,3],
         pch=21, col="black", bg="blue", cex=1.2)
  points(x=inc_e_dec[length(inc_e_dec[,1]),1],
         y=inc_e_dec[length(inc_e_dec[,1]),3],
         pch=21, col="black", bg="red", cex=1.2)

  #sets results file
  unstr_matrix <- S
  unstr_results <- list(data,BdecInc,NBdecInc, new_DI,inc_e_dec,unstr_matrix)
  names(unstr_results) <- c("Original dataset","original TC directions","unstrained TC directions",
                            "unstrained directions and bedding",
                            "inc, E, declination triplets","Unstrain matrix")
  #plot original direction single mode tilt corrected
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(BdecInc,single_mode = TRUE, title="Original directions")

  #plot unstrained directions tilt corrected
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(NBdecInc, single_mode = TRUE, col_d = "red",col_u = "pink", title = "Unstrained directions")

  #save fig
  if(save==TRUE){
    save_pdf(name = paste(name,".pdf"))
    write.csv(unstr_results[[3]],paste(name,"_unstrained TC directions.csv"),
              row.names = FALSE)
    cat("
Unstrained directions saved as csv file.
Figure saved as pdf file.

")
  }
  #restore screen
  par(fig=c(0,1,0,1))
  return(unstr_results)
}

#makes bootstrap stats of EI of initial and final distribution
unstr_stat <- function(unstr_file, nb=1000,confidence=95,hist=TRUE, export=TRUE,name="bootstrap_stat"){
  DI <- unstr_file[[2]]
  NDI <- unstr_file[[3]]
  inc_e_dec <- unstr_file[[5]]
  Inc_E <- as.data.frame(matrix(ncol=3,nrow=0))
  I_E_D0 <- inc_E_finder(DI)
  I_E_Df <- inc_E_finder(NDI)
  NInc_E <- as.data.frame(matrix(ncol=3,nrow=0))

  #bootstrap of original directions
  cat(paste("
",nb, "simulations of original directions:
"))
  for (i in 1:nb) {
    dataprov <- boots_DI(DI)
    I_E_Ed <- inc_E_finder(dataprov)
    I_E_Ed$V1inc <- abs(I_E_Ed$V1inc)
    Inc_E <- rbind(Inc_E, I_E_Ed)
    if(((i%%50)==0)==TRUE) {
      cat(paste(i,"simulations out of",nb,"done
"))
    }
  }
  #bootstrap of unstrained directions
  cat(paste("
",nb, "simulations of unstrained directions:
"))
  for (i in 1:nb) {
    dataprov <- boots_DI(NDI)
    I_E_Ed <- inc_E_finder(dataprov)
    I_E_Ed$V1inc <- abs(I_E_Ed$V1inc)
    NInc_E <- rbind(NInc_E, I_E_Ed)
    if(((i%%50)==0)==TRUE) {
      cat(paste(i,"simulations out of",nb,"done
"))
    }
  }
  #names column
  colnames(Inc_E) <- c("Inc","E","E_dec")
  colnames(NInc_E) <- c("Inc","E","E_dec")
  #makes backup
  Inc_E_bk <- Inc_E
  NInc_E_bk <- NInc_E
  #order for E
  Inc_E <- Inc_E[order(Inc_E$E),]
  NInc_E <- NInc_E[order(NInc_E$E),]
  #trims files for confidence
  conf <- confidence/100
  num <- round((nb*(1-conf))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  Inc_E <- Inc_E[Lconf:Uconf,]
  NInc_E <- NInc_E[Lconf:Uconf,]

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(I_E_D0[1,2]>3.5, ceiling(I_E_D0[1,2]),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #text for figure
  N <- length(DI[,1])
  Inc <- round(inc_e_dec[1,1],digits = 1)
  Ecut <- round(inc_e_dec[1,2],digits = 2)
  V2 <- round(inc_e_dec[1,3],digits = 1)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #unstrained parameters for text2 of figure
  inc_nstr <- round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1)
  E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
  Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

  #put text unstrained in figure
  text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #Plot Bootstrapped data
  points(x=Inc_E$Inc,
         y=Inc_E$E,
         pch=16,
         col=rgb(0.5, 0, 1, 0.15),
         cex=0.6)

  points(x=NInc_E$Inc,
         y=NInc_E$E,
         pch=16,
         col=rgb(1, 0, 0, 0.15),
         cex=0.6)

  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="yellow", lwd=2)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  inc_e_final <- inc_e_dec[length(inc_e_dec[,1]),]
  unx <- inc_e_final[1,1]
  uny <- inc_e_final[1,2]
  points(x=unx,y=uny,
         pch=21, col="black", bg="red", cex=1.2)

  #plot original direction single mode tilt corrected
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(NDI, single_mode = TRUE, col_d = "red",col_u = "pink", title = "Unstrained directions")

  #plot unstrained directions tilt corrected
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(DI,single_mode = TRUE, title="Original directions")

  #plot histogram of Dec if requested
  if(hist==TRUE){
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(Inc_E$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
         breaks = 90, xlab = "", ylab = "",
         main="", cex.axis=0.8, col= "blue", border="blue")
  }

  #unstrained declination
  if(hist==TRUE){
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(NInc_E$E_dec, xlim=c(-90,90), breaks= 90,
         axes=FALSE,xlab="",ylab="",col="red", border="red", main="")
    #plot labels closer than standard to axes
    title(xlab = "Edec(°)", line=1.9, cex=0.3)
    title(ylab = "Frequency", line=2, cex=0.1)
  }
  #plot confidence margin of declination if between -45 and 45
  if(inc_e_final[1,3]>-45 && inc_e_final[1,3]<45){
    Inc_E_Edec <- NInc_E_bk
    Inc_E_Edec <- Inc_E_Edec[order(Inc_E_Edec[,3]),]
    Inc_E_Edec <- Inc_E_Edec[Lconf:Uconf,]
    low_dec <- Inc_E_Edec[1,3]
    up_dec <- Inc_E_Edec[length(Inc_E_Edec[,3]),3]
    abline(v=low_dec,lwd=1,lty=2)
    abline(v=up_dec,lwd=1,lty=2)
  }
  #define low and high inc margins
  inc_conf <- NInc_E_bk
  inc_conf <- inc_conf[order(inc_conf$Inc),]
  inc_conf <- inc_conf[Lconf:Uconf,]
  low_inc <- inc_conf[1,1]
  up_inc <- inc_conf[length(inc_conf[,1]),1]

  if(inc_e_final[1,3]>-45 && inc_e_final[1,3]<45){
    results <- as.data.frame(matrix(ncol= 9, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec","Low_E_dec","High_E_dec")
    results$Inc <- inc_nstr
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- E_nstr
    results$Low_E <- round(min(NInc_E[,2]),digits= 2)
    results$High_E <- round(max(NInc_E[,2]),digits= 2)
    results$E_dec <- Edec_nstr
    results$Low_E_dec <- round(low_dec,digits= 2)
    results$High_E_dec <- round(up_dec,digits=2)
  }else{
    results <- as.data.frame(matrix(ncol= 7, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec")
    results$Inc <- inc_nstr
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- E_nstr
    results$Low_E <- round(min(NInc_E[,2]),digits= 2)
    results$High_E <- round(max(NInc_E[,2]),digits= 2)
    results$E_dec <- Edec_nstr
  }
  print(results, row.names = FALSE)
  par(fig=c(0,1,0,1))
  if(export==TRUE){
    cat("
Results saved as Unstrained_directions_statistic.csv
Graph saved as Unstrained_directions_plot.pdf

")

    #write results file
    write.csv(results,file=paste(name,".csv"),row.names = FALSE)

    #save figure as pdf
    save_pdf(paste(name,".pdf"))
  }
}

# #plot A95 from VGP data and compare with GAPWP
VGP_A95 <- function(VGP,lat=90,long=0,grid=30, auto_cent=TRUE, symbol="c",color="blue",col_A=rgb(1,0,0,0.3), coast=FALSE, on_plot=FALSE, save=FALSE, name="A95",APWP="V23", S_APWP=FALSE){
  library("dplyr", warn.conflicts = FALSE)

  #warning for on-plot, to avoid wrong coordinates
  if(on_plot==TRUE && auto_cent==TRUE) {
    stop("Please SPECIFY center coordinates when on_plot==TRUE",call. = F)
  }

  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}

  #functions spherical (lon=x, lat=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}

  colnames(VGP) <- c("lon","lat")
  vgpsN <- common_DI(VGP,down = ifelse(mean(VGP$lat)<0,FALSE,TRUE))
  PPole <- fisher(vgpsN)
  Plon <- PPole[1,1]
  Plat <- PPole[1,2]
  A <- PPole[1,3]

  #fix point of view
  if(auto_cent==FALSE){
    #center of proj is Lon0 & Lat0
    lon0 <- long
    lat0 <- lat
  }else{
    lon0 <- Plon
    lat0 <- Plat
  }
  #save declination and inc and calculate new system for rotation
  newSlon <- ifelse((Plon+180)>360,Plon-180,Plon+180)
  newSlat <- 90-Plat
  newSlonr <- d2r(newSlon)
  newSlatr <- d2r(newSlat)
  a95 <- A
  circle <- as.data.frame(matrix(ncol=2,nrow=0))
  #loop that create a95 and rotate it around new coordinate (dec, inc)
  for (i in seq(0,360,2)){
    circleP <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(i,(90-a95))
    y <- s2cy(i,(90-a95))
    z <- s2cz(90-a95)
    vec <- as.matrix(c(x,y,z))
    R_elements <- c(cos(newSlatr)*cos(newSlonr), -sin(newSlonr), -sin(newSlatr)*cos(newSlonr),
                    cos(newSlatr)*sin(newSlonr), cos(newSlonr), -sin(newSlatr)*sin(newSlonr),
                    sin(newSlatr), 0, cos(newSlatr))
    R <- matrix(R_elements,nrow=3, byrow=TRUE)
    newvec <- R%*%vec
    newlon <- r2d(atan2(newvec[2,1],newvec[1,1]))
    newlon <- ifelse(newlon<0,newlon+360,newlon)
    #absolute value avoid point outside the graph
    newlat <- r2d(asin(newvec[3,1]))
    circleP[1,1:2] <- c(newlon,newlat)
    circle <- rbind(circle,circleP)
  }
  colnames(circle) <- c("lon","lat")
  circle$x <- c2x(circle$lon,circle$lat)
  circle$y <- c2y(circle$lon,circle$lat)
  circle$cut <- cut(circle$lon,circle$lat)
  #standalone graph or on existing graph
  if (on_plot==FALSE) {sph_ortho(lat = lat0,long = lon0,grid = grid,coast=coast)}

  if(symbol=="c") {sym <- 21}
  else if(symbol=="s") {sym <- 22}
  else if(symbol=="d") {sym <- 23}
  else if(symbol=="t") {sym <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  #plot pole
  Px <- c2x(Plon,Plat)
  Py <- c2y(Plon,Plat)
  Pcut <- cut(Plon,Plat)

  #plot symbol open if behind
  if(Pcut>0) {points(Px,Py,pch=sym,col="black",bg=color)
  }else{points(Px,Py,pch=sym,col="black",bg="white")}

  polygon(circle$x,circle$y, col=col_A, lwd=1.2,lty= ifelse(Pcut>0,1,3))

  if(on_plot==FALSE){
    text <- paste("N: ",round(PPole[1,4],digits=2),"
Long: ", round(PPole[1,1],digits=2),"
Lat: ", round(PPole[1,2],digits=2),"
A95: ", round(PPole[1,3],digits=2))
    text(x=0.75, y=-0.85,pos=4,text, cex= 0.85)
  }

  #plot APWP if requested during process
  pAPWP <- readline("Plot APWP? (y or n): ")
  if(pAPWP=="y"){
    plot_APWP(APWP=APWP, lon0 = lon0,lat0 = lat0,grid = grid,on_plot = T,S_APWP = S_APWP)
  }
  if(save==TRUE){
    save_pdf(name = paste(name,".pdf"),width = 8,height = 8)
    cat("Figure saved as",name, ".pdf")
  }
}

#bootstrap of VGPs
VGP_boot <- function(VGP,nb=1000,lat=90,long=0,grid=30,auto_cent=TRUE,on_plot=FALSE,coast=FALSE,symbol="c",color= "blue",hist=TRUE,text=TRUE,save=FALSE, name="VGP_boot",APWP="V23", S_APWP=FALSE){

  #warning for on-plot, to avoid wrong coordinates
  if(on_plot==TRUE && auto_cent==TRUE) {
    stop("Please SPECIFY center coordinates when on_plot==TRUE",call. = F)
  }

  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #manipulate data
  colnames(VGP) <- c("lon","lat")
  vgpsN <- common_DI(VGP,down = ifelse(mean(VGP$lat)<0,FALSE,TRUE))
  PPole <- fisher(vgpsN)
  Plon <- PPole[1,1]
  Plat <- PPole[1,2]
  #fix point of view
  if(auto_cent==FALSE){
    #center of proj is Lon0 & Lat0
    lon0 <- long
    lat0 <- lat
  }else{
    lon0 <- Plon
    lat0 <- Plat
  }
  if(on_plot==FALSE){sph_ortho(lat=lat0,long=lon0,grid=grid, coast=coast)}
  cat("Bootstrapping.
Simulation ends when", nb, " pseudosamples are saved.

")
  bootlonlat <- as.data.frame(matrix(ncol = 2,nrow = 0))
  n <- 0
  repeat{
    n <- n+1
    VGPb <- boots_DI(VGP)
    VGPb_av <- fisher(VGPb)
    blon <- VGPb_av[1,1]
    blat <- VGPb_av[1,2]
    blonlat <- as.data.frame(t(c(blon,blat)))
    bootlonlat <- rbind(bootlonlat,blonlat)
    x <- c2x(blon,blat)
    y <- c2y(blon,blat)
    cutt <- cut(blon,blat)
    points(x,y,pch=ifelse(cutt>0,16,1),col=rgb(1,0,0,0.15))
    if(((n%%50)==0)==TRUE){

      cat(paste(n,"simulations out of",nb,"done
"))

      if(n==nb) break
    }
  }
  colnames(bootlonlat) <- c("vgp_lon","vgp_lat")
  bootlonlat$Plon <- rep(Plon)
  bootlonlat$Plat <- rep(Plat)
  bootlonlat$delta <- abs(bootlonlat$vgp_lon-bootlonlat$Plon)
  bootlonlat$diff <- r2d(acos((sin(d2r(bootlonlat$vgp_lat))*sin(d2r(bootlonlat$Plat)))+
                                (cos(d2r(bootlonlat$vgp_lat))*cos(d2r(bootlonlat$Plat))*cos(d2r(bootlonlat$delta)))))
  ang_dis <- as.data.frame(bootlonlat$diff)
  ang_dis <- (ang_dis[order(ang_dis[,1]),])
  conf <- 0.95
  Uconf <- round(nb*conf,digits=0)
  angular_conf <- ang_dis[Uconf]

  if(symbol=="c") {sym <- 21}
  else if(symbol=="s") {sym <- 22}
  else if(symbol=="d") {sym <- 23}
  else if(symbol=="t") {sym <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  #plot pole
  Px <- c2x(Plon,Plat)
  Py <- c2y(Plon,Plat)
  Pcut <- cut(Plon,Plat)

  #plot symbol open if behind
  if(Pcut>0) {points(Px,Py,pch=sym,col="black",bg=color)
  }else{points(Px,Py,pch=sym,col="black",bg=NA)}

  #plot angular error estimation
  if(on_plot==FALSE && hist==TRUE){
    par(fig=c(0,0.5,0,0.5), new=TRUE)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    rect(xleft = 0,ybottom = 0,xright = 1,ytop = 1,col = "white",border = "black",cex=0.5)

    par(fig=c(0,0.5,0,0.5), new=TRUE)
    hist(x = bootlonlat$diff,xlab=NA,main="",ylab=NA,
         xlim=c(0,10),breaks=40,cex.axis=0.9,
         col="red",border ="red")
    abline(v=angular_conf,lwd=1,lty=2)
    title(xlab = "Angular distance (°)", line=2, cex=0.2)
    title(ylab = "Frequency", line=2, cex=0.2)
  }
  par(fig=c(0,1,0,1), new=TRUE)
  #plot text with results
  results <- as.data.frame(Plon)
  results$Plat <- Plat
  results$N <- length(VGP$lon)
  results$ang_conf <- angular_conf
  results <- round(results, digits=2)

  if (on_plot==FALSE && text==TRUE){
    text <- paste("N: ",results$N,"
Long: ", results$Plon,"
Lat: ", results$Plat,"
B95: ", results$ang_conf)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

    text(x=0.76, y=0,pos=4,text, cex= 0.85)
  }
  #plot APWP if requested during process
  #plot APWP if requested during process
  pAPWP <- readline("Plot APWP? (y or n): ")
  if(pAPWP=="y"){
    plot_APWP(APWP=APWP, lon0 = lon0,lat0 = lat0,grid = grid,on_plot = T,S_APWP = S_APWP)
  }
  #replot pole data on apwp
  par(fig=c(0,1,0,1), new=TRUE)
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  if(Pcut>0) {points(Px,Py,pch=sym,col="black",bg=color)
  }else{points(Px,Py,pch=sym,col="black",bg=NA)}

  if (save==TRUE){
    write.csv(results, file = paste(name,".csv"), row.names = FALSE)
    save_pdf(name = paste(name,".pdf"),width = 8,height = 8)
    cat("Figure saved as",name, ".pdf
","Result file saved as", name,".csv")
  }
}

#calculate virtual geomagnetic pole(s)
VGP_DI <- function(DI,in_file=FALSE,lat,long,export=TRUE,type="VGPsN",name="VGPs",Prnt=TRUE){
  #conversion functions
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #start data table
  data <- DI[,1:2]
  data <- na.omit(data)
  #add lat and long columns if not present in file
  if(in_file==FALSE){
    data$slat <- rep(lat)
    data$slong <- rep(long)
  }
  #fix col names
  colnames(data) <- c("dec","inc","slat","slong")
  #populate data table with used data
  data$dec_r <- d2r(data$dec)
  data$inc_r <- d2r(data$inc)
  data$slat_r <- d2r(data$slat)
  data$slong_r <- d2r(data$slong)
  #calculate pole colatitude
  data$p_colat_r <- atan(2/tan(data$inc_r))
  data$p_colat_d <- r2d(data$p_colat_r)
  #calculate pole latitude
  data$PLat_r <- asin((sin(data$slat_r)*cos(data$p_colat_r))+
                        (cos(data$slat_r)*sin(data$p_colat_r)*cos(data$dec_r)))
  data$Plat_d<- r2d(data$PLat)
  #Longitudinal difference between site and pole
  data$LDist_r <- asin((sin(data$p_colat_r)*sin(data$dec_r))/cos(data$PLat_r))
  data$LDist_d <- r2d(data$LDist_r)
  #calculate longitude
  data$PLong_d <- ifelse(cos(data$p_colat_r)<(sin(data$slat_r)*sin(data$PLat_r)),
                         data$slong+180-data$LDist_d,
                         data$slong+data$LDist_d)
  #adjust vgp for normal and reversed
  data$Pole_longitude <- ifelse(data$p_colat_d<0,ifelse((data$PLong_d+180)>360,
                                                        data$PLong_d-180,data$PLong_d+180),
                                data$PLong_d)
  data$Pole_latitude <- ifelse(data$p_colat_d<0,-data$Plat_d,data$Plat_d)
  #isolate VGPs with reversals
  VGPs <- data[,16:17]

  #isolate VGPs all normal
  #VGPsN <- data[,c(15,12)]
  VGPsN <- common_DI(VGPs)     #fixed 29.12.2023
  colnames(VGPsN) <- c("Plong_N","Plat_N")

  #calculate average VGP
  PmagPole <- fisher(VGPsN)

  #rename columns of PmagPole
  colnames(PmagPole)[1:2] <- c("long","lat")

  #set coordinates for rotation of pole
  ifelse(PmagPole[1,2]>0,
         RotAZ <- (PmagPole[1,1]+180)%%360,
         RotAZ <-  PmagPole[1,1])
  RotPL <- 90-abs(PmagPole[1,2])

  #rotate VGPs to North pole
  VGPsR <- bed_DI(DI = VGPs,in_file = F,bed_az = RotAZ,
                  bed_plunge = RotPL)
  #if(data[1,3]<0 && PmagPole$lat<0) VGPsR <- flip_DI(VGPsR)
  colnames(VGPsR) <- c("Plong_R","Plat_R")


  if(export==TRUE){
    write.csv(round(PmagPole,digits=2),file=paste(name,"_average_pole.csv"),row.names = F)
    write.csv(round(VGPsN,digits = 2),file=paste(name,"_single_mode.csv"),row.names = F)
    write.csv(round(VGPs, digits = 2),file=paste(name,"_bimodal.csv"),row.names = F)
    write.csv(round(VGPsR,digits = 2),file=paste(name,"_rotated.csv"),row.names = F)
    cat("File exported as csv within the working directory

")
  }
  if(Prnt==TRUE){
    cat("Paleomagnetic pole:

")
    print(round(PmagPole,digits=2),row.names=FALSE)
  }

  if(type=="VGPs"){return(VGPs)}
  if(type=="VGPsN"){return(VGPsN)}
  if(type=="VGPsR"){return(VGPsR)}
}

#plot decl, inc, VGP lat, polarity in stratigraphic depth, and directions and VGP plots if requested
#still under development!!
magstrat_DI <- function(DIP,lat=0,long=0,offset=0,col="red",name="polarity_plot",save=FALSE,plot_ext=TRUE,POLE=TRUE, E.A.=TRUE,cex.main=1,cex.lab=1,cex.axis=1,lwd.grid=1,h_grid=10, Shiny=FALSE){
  library(plyr, warn.conflicts = F)
  library(PmagDiR)
  dat <- na.omit(DIP)
  colnames(dat) <- c("dec","inc","posit")
  #calculate VGPs rotated
  dat[,4:5] <- VGP_DI(dat[,1:2],in_file = F,lat = lat,long = long,export = F,type = "VGPsR",Prnt = F)
  #associated polarity to VGPs, where normal=1, reversed=0
  dat$pol <- ifelse(dat$Plat_R>0,1,0)
  #create reversals empty data frame
  normals <- data.frame(matrix(ncol = 2,nrow = 0))
  colnames(normals) <- c("bottom","top")
  #populate top and bottom of normals table
  for(i in 2:nrow(dat)){
    if((dat[i,6]+dat[i-1,6])==1){
      pos <- (dat[i,3]+dat[i-1,3])/2
      newline <- data.frame(matrix(ncol = 2,nrow = 0))
      ifelse(dat[i,6]==1, newline[1,1] <- pos, newline[1,2] <- pos)
      colnames(newline) <- c("bottom","top")
      normals <- rbind(normals,newline)
    }
  }
  #fill first or last box of normals when empty
  if(is.na(normals[1,1])==TRUE) normals[1,1] <- min(dat$posit)
  if(is.na(normals[nrow(normals),2])==TRUE) normals[nrow(normals),2] <- max(dat$posit)
  #reduce table to lines with top and bottom
  if(nrow(normals)>1){
    for(l in 2:nrow(normals)){
      if(is.na(normals[l-1,2])==TRUE) {normals[l-1,2] <- normals[l,2]}
    }
  }
  #eliminate duplicates
  normals <- na.omit(normals)
  ymin <- plyr::round_any(min(dat$posit), 0.5, f= floor)
  ymax <- plyr::round_any(max(dat$posit), 0.5, f=ceiling)

  #fix declination if offset is not zero
  if(offset!=0){
    offset <- abs(offset)
    dat$dec <- (dat$dec+offset)%%360
    dat$dec <- dat$dec-offset
  }

  ############## PLOT ##############
  #screen splitter matrix
  if(plot_ext==TRUE) {dev.new(width = 10,height = 7,noRStudioGD = T)}
  screen <- matrix(c(1,1,2,2,3,3,4),ncol=7,byrow = T)
  layout(screen)
  inclim <- round(max(abs(dat$inc)), digits = 0)
  #declination
  plot(dat$dec,dat$posit, type="o",
       pch=21,bg=col,ylab="Position (m)",
       xlim=c(-offset,(360-offset)),
       xaxp= c(-offset,(360-offset),4),
       ylim=c(ymin,ymax),
       xlab=NA,
       main="Declination (°)",
       cex.main=cex.main,
       cex.lab=cex.lab,
       cex.axis=cex.axis,
       panel.first= abline(v=c(seq(0-offset,360-offset,90)),
                           h=c(seq(round(min(dat$posit), digits = -1),
                                   round(max(dat$posit), digits = 0),h_grid)),
                           col="gray", lty="dotted",lwd=lwd.grid))
  plot(dat$inc,dat$posit,type="o",
       pch=21,bg=col,ylab=NA,xaxp= c(-90,90,6),
       xlim=c(-90,90),
       ylim=c(ymin,ymax),
       xlab=NA,
       main="Inclination (°)",
       cex.main=cex.main,
       cex.axis=cex.axis,
       panel.first= abline(v=c(seq(-90,90,30)),
                           h=c(seq(round(min(dat$posit), digits = -1),
                                   round(max(dat$posit), digits = 0),h_grid)),
                           col="gray", lty="dotted",lwd=lwd.grid))
  plot(dat$Plat_R,dat$posit,type="o",
       pch=21,bg=col,ylab=NA,xaxp= c(-90,90,4),
       xlim=c(-90,90),
       ylim=c(ymin,ymax),
       xlab=NA,
       main="VGP Lat. (°)",
       cex.main=cex.main,
       cex.axis=cex.axis,
       panel.first= abline(v=c(seq(-90,90,45)),
                           h=c(seq(round(min(dat$posit), digits = -1),
                                   round(max(dat$posit), digits= 0),h_grid)),
                           col="gray", lty="dotted",lwd=lwd.grid))

  #create frame for polarity
  plot(NA,
       xlim=c(0,1), xaxt="n",
       type="n", ylab=NA,
       xlab="",ylim=c(ymin,ymax),
       main="Polarity",
       cex.main=cex.main,
       cex.axis=cex.axis)
  rect(xleft=0,
       ybottom=normals$bottom,
       xright=1,
       ytop=normals$top,
       col=ifelse(nrow(normals==1) && any(dat[,6]==1), "black","white"),
       border=NA)
  if(save==TRUE){
    save_pdf(name =paste(name,".pdf"),width = 10,height = 8)
  }
  if(POLE==TRUE){
    dev.new(width = 7,height = 7,noRStudioGD = T)
    VGPsN <- VGP_DI(dat[,1:2],in_file = F,lat = lat,long = long,export = T,Prnt = T)
    plot_VGP(VGPsN, coast = T, A95 = T,save = save)
  }
  if(E.A.==TRUE){
    dev.new(width = 7,height = 7,noRStudioGD = T)
    plot_DI(dat[,1:2])
    fisher_plot(dat[,1:2],save = save, text=T)
  }
  if(Shiny==TRUE){
    Table_of_normal_polarity_zones <- round(normals,digits=2)
    return(Table_of_normal_polarity_zones)
  }
}


#Temporary IODP plotting script, will be deleted soon
IODP_mag_plot <- function(){
  library("zoo", warn.conflicts = FALSE)
  library("plyr")
  library("dplyr", warn.conflicts = FALSE)

  rm(list = ls()) #clear environment

  ALERT <- readline("IMPORTANT: Works only for data from exp. younger than 362!
Press Enter and select:(1) srmsection file; (2) Core summary file: ")
  if(ALERT == ""){
    srm.data <- read.csv(file.choose())        #import data in csv format
    core.summ <- read.csv(file.choose())       #import core summary file
    core.summ <- na.omit(core.summ)
    srm.data_backup <- srm.data                #copy the file for later use

    #Next two line ask if there is a MTF file to upload
    decornt <- readline("Do you have MTF data?
Type y and select file (downloaded from https://web.iodp.tamu.edu/LORE/) or type n: ")
    drnt.bk <- decornt              #Backup for use in decHist
    if(decornt=="y") {
      orient.file <- read.csv(file.choose())               #if orientation is available, it ask for the file
      loc.dec <- as.numeric(readline("type local declination: "))                  #ask for the local declination
    }

    #Next ask for a discrete declination file if present
    dis.dir <- readline("Do you have a discrete directions file? (y or n)
If yes, load file (.csv in the form depth, dec, inc): ")
    if(dis.dir == "y") {
      dscrt.dirs_all_feat <-  read.csv(file.choose())
      dscrt.dirs <- subset(dscrt.dirs_all_feat,select=c(4,2,3))
      colnames(dscrt.dirs) <- c("depth", "dec", "inc")
    }

    #function generating a dataframe, with one table for each AF step
    new.table <-function(AF){
      depth <-srm.data$Depth.CSF.A..m.[srm.data$Treatment.Value==AF]
      cor <- srm.data$Core[srm.data$Treatment.Value==AF]
      Type <- srm.data$Type[srm.data$Treatment.Value==AF]
      sect <- srm.data$Sect[srm.data$Treatment.Value==AF]
      A.W <- srm.data$A.W[srm.data$Treatment.Value==AF]
      offSet <- srm.data$Offset..cm.[srm.data$Treatment.Value==AF]
      int <- srm.data$Intensity.background...drift.corrected..A.m.[srm.data$Treatment.Value==AF]
      dec <- srm.data$Declination.background...drift.corrected..deg.[srm.data$Treatment.Value==AF]
      inc <- srm.data$Inclination.background...drift.corrected..deg.[srm.data$Treatment.Value==AF]
      new.data <- as.data.frame(cbind(depth, cor, Type, sect, A.W, offSet, int, dec, inc))
    }

    #function eliminating a selected number of points on top and bottom of each core, applied to a data page (e.g. AF0)
    core.TB.filters <- function(x,n.top,n.bot) {
      new.data <- data.frame(depth=numeric(0), cor=numeric(0),
                             int=numeric(0), dec= numeric(0), inc=numeric(0))
      cores.st <- matrix(unique(x$cor))
      for(i in cores.st) {
        data.cores <- filter_all(x, all_vars(x$cor == i))
        data.cores <- data.cores[1:((length(data.cores[,1]))-n.bot),]
        data.cores <- data.cores[n.top:length(data.cores[,1]),]
        new.data <- rbind(new.data,data.cores)
      }
      return(new.data)
    }

    #function copying dec correction values
    MTF.column <- function(x) {
      return(ifelse(x %in% ornt$core, ornt$dec[ornt$core==x], 0.0))
    }

    cores <- matrix(core.summ$Core)                          #List of all cores
    treat <- matrix(sort(unique(srm.data$Treatment.Value)))         #List of all AF steps
    treat.bk <- treat                                               #copy of AF steps list, for later use

    #Next is a loop that stops only when happy about the stratigraphic plots
    repeat {
      srm.data <- srm.data_backup                                         #re-build the original data file (in case some core have been deleted)

      list.of.cores.depth <- subset(core.summ, select=c(4,11,12))
      list.of.cores.depth <- na.omit(list.of.cores.depth)
      colnames(list.of.cores.depth) <- c("Core","Top", "Bottom")
      print(list.of.cores.depth, row.names = FALSE)

      corfilt <- readline(paste("cores range from ", min(cores), " to ", max(cores),     #It state the range of cores and ask if all have to be plotted
                                ". Do you want to plot them all? (y or n): "))
      if(corfilt == "n") {
        upc <- as.numeric(readline("select upper core: "))                 #select a specific cores interval if required, changing the srm.data file
        lowc <- as.numeric(readline("select lower core: "))
        srm.data <- filter_all(srm.data, all_vars(srm.data$Core >= upc ))
        srm.data <- filter_all(srm.data, all_vars(srm.data$Core <= lowc ))
      } else srm.data <- srm.data_backup                                   #rebuild the original srm.data file with all cores

      treat <- treat.bk                                                    #rebuild the original AF steps file in case some has been eliminated during looping
      TreatToPrint <- as.data.frame(treat)
      colnames(TreatToPrint) <- "List of AF Steps"
      print(TreatToPrint, row.names=FALSE)

      #Next allows to eliminate specific AF steps, because in some cases higher field are applied only on limited cores
      AFqst <- readline("Do you want to plot a specific step after NRM? (y or n): ")
      if(AFqst== "y"){
        AFstp <- as.numeric(readline("select step (NRM= 1, second= 2 and so on): "))
      }
      #generate array with a table for any AF steps using the new.table function
      srm.AF.split <- apply(treat, 1, FUN=new.table)

      ##########Next generate the parameters for the log figures, only for the NRM and the last AF step

      AF0 <- as.data.frame(srm.AF.split[[1]])                   #Table NRM
      AF0[,1] <- as.numeric(AF0[,1])                            #convert chr in num
      AF0[,2] <- as.numeric(AF0[,2])
      AF0[,4] <- as.numeric(AF0[,4])
      AF0[,6] <- as.numeric(AF0[,6])
      AF0[,7] <- as.numeric(AF0[,7])
      AF0[,8] <- as.numeric(AF0[,8])
      AF0[,9] <- as.numeric(AF0[,9])
      AF0 <- AF0[order(AF0$depth),]

      last.AF <- ifelse(AFqst== "y",treat[AFstp], treat[length(treat)])    #last AF step
      last.page <- ifelse(AFqst== "y", AFstp, length(treat))                                #index of the last AF step
      ymin <- round_any(min(AF0$depth), 10, f= floor)           #min depth of columns, approximated by 10 meters
      ymax <- round_any(max(AF0$depth), 10, f=ceiling)          #max depth of columns, approximated by 10 meters
      ysubs <- ((ymax-ymin))/10                                 #Subdivision of depth scale in 10 meters
      AF.last <- as.data.frame(srm.AF.split[[last.page]])       #Table last AF step
      AF.last[,1] <- as.numeric(AF.last[,1])                    #convert chr in num
      AF.last[,2] <- as.numeric(AF.last[,2])
      AF.last[,4] <- as.numeric(AF.last[,4])
      AF.last[,6] <- as.numeric(AF.last[,6])
      AF.last[,7] <- as.numeric(AF.last[,7])
      AF.last[,8] <- as.numeric(AF.last[,8])
      AF.last[,9] <- as.numeric(AF.last[,9])
      AF.last <- AF.last[order(AF.last$depth),]

      #Next eliminates a number of selected point measurement either on top or bottom of each core, if required
      core.TB <- readline("filtering top and bottom of single cores? (y or n): ")

      if(core.TB =="y"){
        n.top <- 1+as.numeric(readline("how many points on top?  "))
        n.bot <- as.numeric(readline("how many point at bottom? "))
        AF0 <- core.TB.filters(AF0,n.top = n.top, n.bot=n.bot)
        AF.last <- core.TB.filters(AF.last,n.top = n.top, n.bot=n.bot)
      }

      NRM <- log10(AF0$int)                                     #NRM intensity on logarithmic scale
      NRM.last <- log10(AF.last$int)                            #NRM after demag
      dec.NRM <- as.numeric(AF0$dec)                            #NRM declination
      dec.last <- as.numeric(AF.last$dec)                       #declination after last AF step

      dscrt.qst <- ifelse(dis.dir=="y", readline("do you want to plot discrete samples directions? (y or n): "), "n")

      decqst <- readline("Do you want to plot Declination? (y or n): ")         #if declination plot is required, it split the figure screen in five, otherwise 3

      m3 <- matrix(c(1,1,2,2,2,3,3,3,4,4,4,5,5), ncol=13, byrow=TRUE)

      m5 <- matrix(c(1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7), ncol= 19, byrow= TRUE)

      #open new device for plotting figure
      dev.new(width = 11,height = 9,noRStudioGD = T)

      repeat {                                      #repeat the whole cycle because if the running mean of inclination is not ok, it can be redrawn
        ifelse(decqst == "y", layout(m5), layout(m3))               #split screen

        ############# create column with cores name and depth
        cornames.4.plot <- list.of.cores.depth
        cornames.4.plot$mean <- (cornames.4.plot[,2]+cornames.4.plot[,3])/2
        cornames.4.plot$xpos <- 0.5

        plot(0,                         #create frame
             ylim=c(ymax, ymin),
             xlim=c(0,1),
             yaxp= c(ymax, ymin, ysubs),
             xaxt="n",
             type="n",
             main="Core",
             ylab="Depth (m CSF-A)",
             xlab="")
        rect(xleft=0,
             ybottom=cornames.4.plot$Bottom,
             xright=1,
             ytop=cornames.4.plot$Top,
             col="gray")
        text(x=cornames.4.plot$xpos,
             y=cornames.4.plot$mean,
             labels=cornames.4.plot$Core)
        ###############################

        plot(x= NRM, y= AF0$depth, type= "p",                     #plot NRM and last AF intensity
             pch = 16,
             col="blue",
             cex= 0.7,
             xlim=c(min(log10(srm.data$Intensity.background...drift.corrected..A.m.)),
                    max(log10(srm.data$Intensity.background...drift.corrected..A.m.))),
             ylim=c(ymax, ymin),
             main = "NRM
intensity (log A/m)",
             xlab= "",
             ylab= "",
             yaxp= c(ymax, ymin, ysubs))
        points(x=NRM.last,
               y= AF.last$depth,
               pch = 16,
               col="cyan",
               cex= 0.7)

        if(decqst == "y"){                                       #if declination plot is required, it ask for a orientation file

          #################  all that follows add two columns (5, 6) to the pages with the NRM (AF0) and last AF step (AF.last), 5= correction, 6= final dec.
          #################  in 6 the magnetometer dec is copied if there is no correction applicable

          if(decornt == "n") {       #plot NRM declination if no orientation is given
            plot(x= AF0$dec,
                 y= AF0$depth,
                 type= "p",
                 pch = 16,
                 col="blue",
                 cex= 0.7,
                 ylim=c(ymax, ymin),
                 main= "NRM
declination (°)",
                 xlab= "",
                 ylab= "",
                 xaxp= c(0, 360, 4),
                 yaxp= c(ymax, ymin, ysubs))

          } else {
            ornt <- subset(orient.file, select=c(4,14))      #Isolate the declination correction
            colnames(ornt) <- c("core", "dec")


            cor.num_AF0 <- as.matrix(AF0$cor)                     #column with core number
            cor.num_AF.last <- as.matrix(AF.last$cor)

            new_column_MTF_AF0 <- data.frame(apply(cor.num_AF0, MARGIN = 1, MTF.column))  #generate the column with the dec correction
            colnames(new_column_MTF_AF0) <- "MTF"
            new_column_MTF_AF0$IGRF <- loc.dec
            new_column_MTF_AF0$MTF.IGRF <- new_column_MTF_AF0[,1]+new_column_MTF_AF0[,2]
            AF0 <- cbind(AF0,new_column_MTF_AF0[,3])
            colnames(AF0) <- c("depth", "cor","Type","sect","A.W", "offSet","int", "dec","inc","MTF+IGRF")
            AF0$final.dec <- ifelse((AF0$dec+AF0$`MTF+IGRF`) >= 360,
                                    (AF0$dec+AF0$`MTF+IGRF`) - 360,
                                    (AF0$dec+AF0$`MTF+IGRF`))
            AF0$final.dec <- ifelse(AF0$final.dec<0, AF0$final.dec+360, AF0$final.dec)

            new_column_MTF_AF.last <- data.frame(apply(cor.num_AF.last, MARGIN = 1, MTF.column))
            colnames(new_column_MTF_AF.last) <- "MTF"
            new_column_MTF_AF.last$IGRF <- loc.dec
            new_column_MTF_AF.last$MTF.IGRF <- new_column_MTF_AF.last[,1]+new_column_MTF_AF.last[,2]
            AF.last <- cbind(AF.last,new_column_MTF_AF.last[,3])
            colnames(AF.last) <- c("depth", "cor","Type","sect","A.W", "offSet","int", "dec","inc","MTF+IGRF")
            AF.last$final.dec <- ifelse((AF.last$dec+AF.last$`MTF+IGRF`) >= 360,
                                        (AF.last$dec+AF.last$`MTF+IGRF`) - 360,
                                        (AF.last$dec+AF.last$`MTF+IGRF`))
            AF.last$final.dec <- ifelse(AF.last$final.dec<0, AF.last$final.dec+360,AF.last$final.dec)


            ########## end of calculation and columns compilation ##########################

            plot(x= AF0$final.dec,       #Plot NRM corrected declination
                 y= AF0$depth,
                 type= "p",
                 pch = 16,
                 col="blue",
                 cex= 0.7,
                 ylim=c(ymax, ymin),
                 main= "NRM
declination (°)",
                 xlab= "",
                 ylab= "",
                 xaxp= c(0, 360, 4),
                 yaxp= c(ymax, ymin, ysubs))
          }
        }

        plot(x= AF0$inc,       #plot NRM inclination
             y= AF0$depth,
             type= "p",
             pch = 16,
             col="blue",
             cex= 0.7,
             ylim=c(ymax, ymin),
             xlim=c(-90,90),
             main= "NRM
inclination (°)",
             xlab= "",
             ylab= "",
             xaxp= c(-90, 90, 4),
             yaxp= c(ymax, ymin, ysubs))
        abline(v=0,lwd=1, lty=2)

        if(decqst == "y"){                               #condition when declination plot is requested or not

          if(decornt == "n") {
            plot(x= AF.last$dec, y= AF.last$depth, type= "p",          #plot dec of last AF without orientation
                 pch = 16,
                 col="cyan",
                 cex= 0.7,
                 ylim=c(ymax, ymin),
                 main= paste(last.AF,"mT
","Declination (°)", sep = ""),
                 xlab= "",
                 ylab= "",
                 xaxp= c(0, 360, 4),
                 yaxp= c(ymax, ymin, ysubs))

            if(dscrt.qst == "y") points(x=dscrt.dirs$dec, y=dscrt.dirs$depth,     #ask and plot discrete declination
                                        type= "p",
                                        pch=21,
                                        col="black",
                                        bg="red")

          } else {
            plot(x= AF.last$final.dec, y= AF.last$depth, type= "p",          #Plot declination with correction
                 pch = 16,
                 col="cyan",
                 cex= 0.7,
                 ylim=c(ymax, ymin),
                 main= paste(last.AF,"mT
","Declination (°)", sep = ""),
                 xlab= "",
                 ylab= "",
                 xaxp= c(0, 360, 4),
                 yaxp= c(ymax, ymin, ysubs))

            if(dscrt.qst == "y") points(x=dscrt.dirs$dec, y=dscrt.dirs$depth,    #ask and plot discrete declination
                                        type= "p",
                                        pch=21,
                                        col="black",
                                        bg="red")
          }
        }
        plot(x= AF.last$inc,                 #Plot last AF inclination
             y= AF.last$depth, type= "p",
             pch = 16,
             col= "cyan",
             cex= 0.7,
             ylim=c(ymax, ymin),
             xlim=c(-90,90),
             main= paste(last.AF,"mT
","Inclination (°)", sep = ""),
             xlab= "",
             ylab= "",
             xaxp= c(-90, 90, 4),
             yaxp= c(ymax, ymin, ysubs))
        abline(v=0,lwd=1, lty=2)

        if(dscrt.qst == "y") points(x=dscrt.dirs$inc, y=dscrt.dirs$depth,    #it adds the discrete samples inclination
                                    type= "p",
                                    pch=21,
                                    col="black",
                                    bg="red")


        ########## Next give the opportunity to add a running mean to inclination data

        rnmean <- readline("plotting also running mean of inclination? (y or n): ")

        if (rnmean== "y") {                          #if running mean is requested, it ask for the number of points to average
          pt <- as.numeric(readline("number of point you want to use for the running mean: "))

          depth.Inc.last <- subset(AF.last, select=c(1,9))
          depth.Inc.last <- depth.Inc.last[order(depth.Inc.last$depth),]
          aver.CSF.A <- rollmean(as.numeric(depth.Inc.last$depth), k= pt)       #moving average of depth
          aver.inc <- rollmean(as.numeric(depth.Inc.last$inc), k= pt)         #moving average of inc

          points(x= aver.inc, y= aver.CSF.A, type= "l", lwd=1.5, col="black")               #it adds the running mean to the log
          runmeanqst <- readline("Do you want to change number of averaged points? (y or n): ")
          if (runmeanqst != "y") break
        }
        if (rnmean != "y") break
      }

      ######### Next plot a polarity indication #####

      DecOrInc <- readline("Do you want to use Dec or Inc for interpreting polarity? (d or i): ")
      if (DecOrInc == "i") {
        N.S.EM <- readline("Northern or southern emisphere? (n or s): ")
        if (rnmean == "y") {
          pol.run.mean <- readline("Do you want to use running mean of inclination? (y or n): ")
          if (pol.run.mean =="n") {
            if (N.S.EM == "n") AF.last$N <- ifelse(AF.last$inc>0,1,0)
            if (N.S.EM == "s") AF.last$N <- ifelse(AF.last$inc>0,0,1)
            only.N <- AF.last[,c(1,10)]
          }
          if (pol.run.mean =="y"){
            if (N.S.EM == "n") roll.mean.pol <- ifelse(aver.inc>0,1,0)
            if (N.S.EM == "s") roll.mean.pol <- ifelse(aver.inc>0,0,1)
            roll.mean.plot <- as.data.frame(cbind(aver.CSF.A, roll.mean.pol))
            colnames(roll.mean.plot) <- c("depth", "N")
            only.N <- roll.mean.plot
          }
        }
        if (rnmean=="n") {
          if (N.S.EM == "n") AF.last$N <- ifelse(AF.last$inc>0,1,0)
          if (N.S.EM == "s") AF.last$N <- ifelse(AF.last$inc>0,0,1)
          only.N <- AF.last[,c(1,10)]
        }
      }
      if (DecOrInc == "d") {
        l.dec = as.numeric(readline("Set declination reversed lower angle: "))
        u.dec = as.numeric(readline("Set declination reversed higher angle: "))
        if (decornt=="n"){AF.last$final.dec <- AF.last$dec}
        AF.last$N <- ifelse(AF.last$final.dec>l.dec, ifelse(AF.last$final.dec<u.dec,0,1),1)
        only.N <- AF.last[,c(1,10)]
      }

      #create reversals empty data frame
      normals <- data.frame(matrix(ncol = 2,nrow = 0))
      colnames(normals) <- c("bottom","top")
      #populate top and bottom of normals table
      for(i in 2:nrow(only.N)){
        if((only.N[i,2]+only.N[i-1,2])==1){
          pos <- (only.N[i,1]+only.N[i-1,1])/2
          newline <- data.frame(matrix(ncol = 2,nrow = 0))
          ifelse(only.N[i,2]==1, newline[1,1] <- pos, newline[1,2] <- pos)
          colnames(newline) <- c("bottom","top")
          normals <- rbind(normals,newline)
        }
      }
      #fill first or last box of normals when empty
      if(is.na(normals[1,1])==TRUE) normals[1,1] <- min(only.N[,1])
      if(is.na(normals[nrow(normals),2])==TRUE) normals[nrow(normals),2] <- max(only.N[,1])
      #reduce table to lines with top and bottom
      if(nrow(normals)>1){
        for(l in 2:nrow(normals)){
          if(is.na(normals[l-1,2])==TRUE) {normals[l-1,2] <- normals[l,2]}
        }
      }
      #eliminate duplicates


      plot(0,                         #create frame
           ylim=c(ymax, ymin),
           xlim=c(0,1),
           yaxp= c(ymax, ymin, ysubs),
           xaxt="n",
           type="n",
           main="polarity",
           ylab="",
           xlab="")
      rect(xleft=0,
           ybottom=normals$bottom,
           xright=1,
           ytop=normals$top,
           col=ifelse(nrow(normals==1) && any(only.N[,2]==1), "black","white"),
           border=NA)

      # arrows(x0=0,
      #        y0=only.N$depth,
      #        x1=0+(only.N$N),
      #        y1=only.N$depth,
      #        code=3,
      #        length=0,
      #        lwd=0.5,
      #        col="black")

      redo <- readline("are you happy? (y or n): ")
      if (redo=="y") rm(only.N)
      if(redo=="y") break                                                  #break the loop if happy, or it start from beginning
    }

    rm(decornt,dis.dir)                                                    #delete initial questions, otherwise they affect next elaborations

    write.csv(AF.last, file="Table_data.csv", row.names = FALSE)

    qst <- readline("want histogram of inc? (y or n): ")                   #ask for histogram of inclination

    if(qst=="y") {
      repeat {                                                             #loop for hist, it breaks if happy
        par(mfrow=c(1,1))
        if(qst=="y") subs <- as.numeric(readline("insert the bin size in degrees: "))
        bin <- 90/subs
        hist(AF.last$inc, plot=TRUE,
             xlim= c(-90, 90),
             xaxp= c(-90, 90, 6),
             xlab= "Inclination (°)",
             ylab= paste("N. of directions (Total= ",length(AF.last$inc), ")"),
             main=paste(max(treat),"mT
","Inclination", sep = ""),
             col= "red",
             breaks= bin)
        redoh <-readline("are you happy? (y or n): ")
        if(redoh=="y") break                                              #breaks the loop if happy
      }

    }
    qstd <- readline("want histogram of Dec? (y or n): ")
    if(qstd=="y") {
      repeat {                                                             #loop for hist, it breaks if happy
        AF.last.bk <- AF.last
        par(mfrow=c(1,1))
        subs <- as.numeric(readline("insert the bin size in degrees: "))
        bin <- 360/subs
        if (drnt.bk=="n"){AF.last.bk$final.dec <- AF.last.bk$dec}
        final_dec <- AF.last.bk$final.dec
        final_dec <- ifelse(final_dec>270, final_dec-360,final_dec)
        hist(final_dec, plot=TRUE,
             xlim= c(-90, 270),
             xaxp= c(-90, 270, 4),
             xlab= "Declination (°)",
             ylab= paste("N. of directions (Total= ",length(AF.last.bk$final.dec), ")"),
             main=paste(max(treat),"mT
","Declination", sep = ""),
             col= "blue",
             breaks= bin)
        redohd <-readline("are you happy? (y or n): ")
        if(redohd=="y") break                                              #breaks the loop if happy
      }
    }
  }
}




art_make( BASENAME_ONLY MODULE_LIBRARIES
  lbne-raw-data_Overlays
  dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
  lbne_daqinput35t       #For online to offline channel map
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  dune_DuneObj
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
  lardata_ArtDataHelper
  nusimdata_SimulationBase	
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}

  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  EXCLUDE LinHitFinderAlg1_tool.cc LinHitFinderAlg2_tool.cc
  )

simple_plugin(LinHitFinderAlg1 "tool"
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)

simple_plugin(LinHitFinderAlg2 "tool"
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)

add_subdirectory(HitDumper)

install_fhicl()
install_headers()
install_source()

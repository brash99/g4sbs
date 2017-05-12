# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# compile C with /usr/bin/cc
C_FLAGS = -O3 -DNDEBUG   -DCTEQ_TBL_PATH=\"/home/meriem/g4sbs_detecRTPC/g4sbs/cteq-tbls\"

C_DEFINES = -DG4INTY_USE_QT -DG4INTY_USE_XT -DG4UI_USE -DG4UI_USE_QT -DG4UI_USE_TCSH -DG4VERBOSE -DG4VIS_USE -DG4VIS_USE_OPENGL -DG4VIS_USE_OPENGLQT -DG4VIS_USE_OPENGLX -DG4_STORE_TRAJECTORY

C_INCLUDES = -I/home/meriem/g4sbs_detecRTPC/g4sbs/include -I/home/meriem/root/include -I/home/meriem/g4sbs_detecRTPC/g4sbs/src -I/home/meriem/g4sbs_detecRTPC/g4sbs/src/dss2007 -I/home/meriem/g4sbs_detecRTPC/g4sbs/build/include -isystem /home/meriem/geant4.10.03.p01/source/analysis/g4tools/include -isystem /home/meriem/geant4.10.03.p01/source/analysis/accumulables/include -isystem /home/meriem/geant4.10.03.p01/source/analysis/csv/include -isystem /home/meriem/geant4.10.03.p01/source/analysis/hntools/include -isystem /home/meriem/geant4.10.03.p01/source/analysis/management/include -isystem /home/meriem/geant4.10.03.p01/source/analysis/root/include -isystem /home/meriem/geant4.10.03.p01/source/analysis/xml/include -isystem /home/meriem/geant4.10.03.p01/source/digits_hits/detector/include -isystem /home/meriem/geant4.10.03.p01/source/digits_hits/digits/include -isystem /home/meriem/geant4.10.03.p01/source/digits_hits/hits/include -isystem /home/meriem/geant4.10.03.p01/source/digits_hits/scorer/include -isystem /home/meriem/geant4.10.03.p01/source/digits_hits/utils/include -isystem /home/meriem/geant4.10.03.p01/source/error_propagation/include -isystem /home/meriem/geant4.10.03.p01/source/event/include -isystem /home/meriem/geant4.10.03.p01/source/externals/clhep/include -isystem /home/meriem/geant4.10.03.p01/source/externals/zlib/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/biasing/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/divisions/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/magneticfield/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/management/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/navigation/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/solids/Boolean/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/solids/CSG/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/solids/specific/include -isystem /home/meriem/geant4.10.03.p01/source/geometry/volumes/include -isystem /home/meriem/geant4.10.03.p01/source/global/HEPGeometry/include -isystem /home/meriem/geant4.10.03.p01/source/global/HEPNumerics/include -isystem /home/meriem/geant4.10.03.p01/source/global/HEPRandom/include -isystem /home/meriem/geant4.10.03.p01/source/global/management/include -isystem /home/meriem/geant4.10.03.p01/source/graphics_reps/include -isystem /home/meriem/geant4.10.03.p01/source/intercoms/include -isystem /home/meriem/geant4.10.03.p01/source/interfaces/GAG/include -isystem /home/meriem/geant4.10.03.p01/source/interfaces/basic/include -isystem /home/meriem/geant4.10.03.p01/source/interfaces/common/include -isystem /home/meriem/geant4.10.03.p01/source/materials/include -isystem /home/meriem/geant4.10.03.p01/source/parameterisations/gflash/include -isystem /home/meriem/geant4.10.03.p01/source/particles/adjoint/include -isystem /home/meriem/geant4.10.03.p01/source/particles/bosons/include -isystem /home/meriem/geant4.10.03.p01/source/particles/hadrons/barions/include -isystem /home/meriem/geant4.10.03.p01/source/particles/hadrons/ions/include -isystem /home/meriem/geant4.10.03.p01/source/particles/hadrons/mesons/include -isystem /home/meriem/geant4.10.03.p01/source/particles/leptons/include -isystem /home/meriem/geant4.10.03.p01/source/particles/management/include -isystem /home/meriem/geant4.10.03.p01/source/particles/shortlived/include -isystem /home/meriem/geant4.10.03.p01/source/particles/utils/include -isystem /home/meriem/geant4.10.03.p01/source/persistency/ascii/include -isystem /home/meriem/geant4.10.03.p01/source/persistency/mctruth/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/builders/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/decay/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/electromagnetic/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/factory/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/gamma_lepto_nuclear/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/hadron_elastic/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/hadron_inelastic/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/ions/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/limiters/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/constructors/stopping/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/lists/include -isystem /home/meriem/geant4.10.03.p01/source/physics_lists/util/include -isystem /home/meriem/geant4.10.03.p01/source/processes/biasing/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/biasing/generic/include -isystem /home/meriem/geant4.10.03.p01/source/processes/biasing/importance/include -isystem /home/meriem/geant4.10.03.p01/source/processes/cuts/include -isystem /home/meriem/geant4.10.03.p01/source/processes/decay/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/adjoint/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/dna/processes/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/dna/models/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/dna/utils/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/dna/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/dna/molecules/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/dna/molecules/types/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/highenergy/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/lowenergy/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/muons/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/pii/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/polarisation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/standard/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/utils/include -isystem /home/meriem/geant4.10.03.p01/source/processes/electromagnetic/xrays/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/cross_sections/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/abla/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/abrasion/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/binary_cascade/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/cascade/cascade/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/coherent_elastic/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/ablation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/evaporation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/fermi_breakup/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/fission/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/gem_evaporation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/handler/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/multifragmentation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/photon_evaporation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/de_excitation/util/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/em_dissociation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/fission/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/im_r_matrix/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/inclxx/utils/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/inclxx/incl_physics/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/inclxx/interface/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/lend/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/lepto_nuclear/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/particle_hp/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/parton_string/diffraction/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/parton_string/hadronization/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/parton_string/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/parton_string/qgsm/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/pre_equilibrium/exciton_model/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/qmd/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/radioactive_decay/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/quasi_elastic/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/rpg/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/theo_high_energy/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/models/util/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/processes/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/stopping/include -isystem /home/meriem/geant4.10.03.p01/source/processes/hadronic/util/include -isystem /home/meriem/geant4.10.03.p01/source/processes/management/include -isystem /home/meriem/geant4.10.03.p01/source/processes/optical/include -isystem /home/meriem/geant4.10.03.p01/source/processes/phonon/include -isystem /home/meriem/geant4.10.03.p01/source/processes/parameterisation/include -isystem /home/meriem/geant4.10.03.p01/source/processes/scoring/include -isystem /home/meriem/geant4.10.03.p01/source/processes/transportation/include -isystem /home/meriem/geant4.10.03.p01/source/readout/include -isystem /home/meriem/geant4.10.03.p01/source/run/include -isystem /home/meriem/geant4.10.03.p01/source/track/include -isystem /home/meriem/geant4.10.03.p01/source/tracking/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/FukuiRenderer/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/HepRep/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/RayTracer/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/Tree/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/VRML/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/XXX/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/externals/gl2ps/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/gMocren/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/management/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/modeling/include -isystem /home/meriem/geant4.10.03.p01/source/visualization/OpenGL/include 


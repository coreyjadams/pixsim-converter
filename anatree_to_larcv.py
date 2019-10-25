import numpy
from larcv import larcv
import glob
import pickle
import os

import root_numpy

import argparse

angle = 35.7
# angle = 89
cos_theta = numpy.cos(numpy.deg2rad(angle))
sin_theta = numpy.sin(numpy.deg2rad(angle))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', help="Name of file")
    args = parser.parse_args()


    input_anatree_file_name = args.filename

    # convert_training_set('')
    convert_file(input_anatree_file_name)

# files = files[0:5]

def get_dune_meta3D():

    dune_meta = larcv.ImageMeta3D()
    # set_dimension(size_t axis, double image_size, size_t number_of_voxels, double origin = 0);

    pixel_size = 0.4
    dune_meta.set_dimension(0, 360, int(360./pixel_size), 0)
    dune_meta.set_dimension(1, 200, int(200./pixel_size), -100)
    dune_meta.set_dimension(2, 500, int(500./pixel_size),  0)


    return dune_meta

def get_dune_meta2D(plane):

    dune_meta = larcv.ImageMeta2D()
    dune_meta.set_projection_id(plane)
    # set_dimension(size_t axis, double image_size, size_t number_of_voxels, double origin = 0);

    pixel_size = 0.4

    if plane == 2:
        image_size = 500
        origin = 0
    if plane == 1:
        xmin = cos_theta*0   - sin_theta*100.
        xmax = cos_theta*500 - sin_theta*-100
        # print("xmin: ", xmin)
        # print("xmax: ", xmax)
        # cos_theta * ymin
        image_size = xmax - xmin
        origin = xmin
    if plane == 0:
        xmin = cos_theta*0    + sin_theta*-100.
        xmax = cos_theta*500  + sin_theta*100
        # print("xmin: ", xmin)
        # print("xmax: ", xmax)
        image_size = xmax - xmin
        origin = xmin

    number_of_voxels = int(image_size / pixel_size)
    # print("number_of_voxels: ", number_of_voxels)
    # print("image_size: ", image_size)

    if image_size != int(image_size):
        # Add an extra voxel:
        number_of_voxels = int( (image_size / pixel_size) + 1)

        #Recompute the image size to be exact:
        image_size = number_of_voxels*pixel_size
        # print("image_size: ", image_size)
        number_of_voxels = image_size / pixel_size
        # print("number_of_voxels: ", number_of_voxels)
        number_of_voxels = int(image_size / pixel_size)
        # print("number_of_voxels: ", number_of_voxels)



    # if plane == 0:
    #     pixel_size = 0.467
    # else:
    #     pixel_size = 0.479
    dune_meta.set_dimension(0,  image_size, number_of_voxels, origin)
    dune_meta.set_dimension(1,  360, int(360./pixel_size), 0)

    print(dune_meta.dump())

    return dune_meta


def convert_file(input_anatree_file_name):

    # First, open the anatree file with root_numpy and convert it:
    events = root_numpy.root2array(input_anatree_file_name,"anatree/anatree")


    # Second, initialize larcv:
    io_manager = larcv.IOManager(larcv.IOManager.kWRITE)

    output = input_anatree_file_name.replace(".root", ".h5")
    io_manager.set_out_file(output)
    io_manager.initialize()


    # We do several steps: converting particle information,
    # Then convert IDEs into cluster3D, then project them into 2D on 3 planes.

    dune_meta_3D = get_dune_meta3D()
    dune_meta_2D = [ get_dune_meta2D(i) for i in [0, 1, 2]]

    for i_event in range(len(events)):
        print("Converting entry ", i_event)
        event = events[i_event]
        run    = event['run']
        subrun = event['subrun']
        evt    = event['event']

        io_manager.set_id(int(run), int(subrun), int(evt))

        # Store the neutrino information:
        larcv_neut_particle = larcv.EventParticle.to_particle(io_manager.get_data("particle","neutrino"))
        particle = larcv.Particle()

        particle.pdg_code(int(event['nu_pdg'][0]))
        # particle.(event['nu_ndau'][0])
        particle.nu_current_type(int(event['nu_ccnc'][0]))
        particle.nu_interaction_type(int(event['nu_mode'][0]))
        # particle.(int(event['nu_inttype'][0]))
        particle.position(
            event['nu_StartPointx'][0],
            event['nu_StartPointy'][0],
            event['nu_StartPointz'][0], 0.0
        )
        particle.end_position(
            event['nu_EndPointx'][0],
            event['nu_EndPointy'][0],
            event['nu_EndPointz'][0], 0.0
        )
        particle.first_step(
            event['nu_Vertexx'][0],
            event['nu_Vertexy'][0],
            event['nu_Vertexz'][0], 0.0
        )
        particle.momentum(
            event['nu_StartPx'][0],
            event['nu_StartPy'][0],
            event['nu_StartPz'][0]
        )
        larcv_neut_particle.emplace_back(particle)
        # particle.(event['nu_EndPx'][0])
        # particle.(event['nu_EndPy'][0])
        # particle.(event['nu_EndPz'][0])

        # Get the geant particle information:
        n_particles = event['geant_list_size']
        larcv_particle = larcv.EventParticle.to_particle(io_manager.get_data("particle","segment"))
        for p in range(n_particles):
            particle = larcv.Particle()

            particle.pdg_code(int(event['PDG'][p]))
            particle.energy_init(event['StartEnergy'][p])
            particle.position(
                event['StartPointx'][p],
                event['StartPointy'][p],
                event['StartPointz'][p], 0.0
            )
            particle.end_position(
                event['EndPointx'][p],
                event['EndPointy'][p],
                event['EndPointz'][p], 0.0
            )
            particle.momentum(
                event['StartPx'][p],
                event['StartPy'][p],
                event['StartPz'][p]
            )
            particle.creation_process(str(event['Process'][p]))

            particle.track_id(int(event['TrackId'][p]))
            particle.parent_track_id(int(event['Mother'][p]))
            particle.ancestor_track_id(int(event['isPrimary'][p]))

            # particle.(event['nu_EndPx'][p])
            # particle.(event['nu_EndPy'][p])
            # particle.(event['nu_EndPz'][p])
            larcv_particle.emplace_back(particle)


        # Get the IDE information:
        ides_x  = event['ides_x']
        ides_y  = event['ides_y']
        ides_z  = event['ides_z']
        ides_e  = event['ides_energy']
        ides_id = event['ides_tid']
        ides_ne = event['ides_numElectrons']

        # print()
        # print("x min:",  numpy.min(ides_x))
        # print("y min:",  numpy.min(ides_y))
        # print("z min:",  numpy.min(ides_z))
        # print("x max:",  numpy.max(ides_x))
        # print("y max:",  numpy.max(ides_y))
        # print("z max:",  numpy.max(ides_z))
        # # print(ides_pixel_y[50:])
        # # exit()

        # print(numpy.unique(ides_x))

        # Create a sparse cluster object and populate it with the ides:
        ev_clusters3d = larcv.EventSparseCluster3D.to_sparse_cluster(io_manager.get_data("cluster3d", "segment"))
        ev_clusters2d = larcv.EventSparseCluster2D.to_sparse_cluster(io_manager.get_data("cluster2d", "segment"))


        clusters_3d = larcv.SparseCluster3D()
        clusters_2d = [ larcv.SparseCluster2D() for i in [0, 1, 2] ]
        # Create a cluster for each of the particles in the loop above, plus one for anything that
        # can't be found.
        clusters_3d.meta(dune_meta_3D)
        _ = [clusters_2d[i].meta(dune_meta_2D[i]) for i in [0,1,2]]


        track_id_map = {}
        cluster_vec_3d = larcv.VectorOfVoxelSet()
        cluster_vec_2d = [ larcv.VectorOfVoxelSet() for i in [0,1,2]]

        cluster_vec_3d.resize(len(event['TrackId']) + 1)
        _ = [cluster_vec_2d[i].resize(len(event['TrackId']) + 1 ) for i in [0,1,2]]

        for i in range(len(event['TrackId'])):
            track_id_map[event['TrackId'][i]] = i
        unknown_index = len(event['TrackId'])

        coords_3d = larcv.VectorOfDouble()
        coords_2d = larcv.VectorOfDouble()
        coords_3d.resize(3)
        coords_2d.resize(2)
        for i in range(len(ides_x)):

            coords_3d[0] = ides_x[i]
            coords_3d[1] = ides_y[i]
            coords_3d[2] = ides_z[i]
            index = dune_meta_3D.position_to_index(coords_3d)

            # for _temp in [0,1,2]:
            #     position = dune_meta_3D.position(index, _temp)
            #     coordinate = dune_meta_3D.coordinate(index,_temp)
            #     print()
            #     print(coords_3d[_temp], " vs ", position)
            #     if _temp == 1:
            #         print(coords_3d[_temp], coords_3d[_temp] / 0.4, dune_meta_3D.position_to_coordinate(coords_3d[_temp], _temp))
            #         # print("Y: ", position, " - ", ides_pixel_y[i], " = ", position - ides_pixel_y[i], coords_3d[_temp] / 0.4, coordinate)
            #     elif _temp == 2:
            #         print(coords_3d[_temp], coords_3d[_temp] / 0.4, dune_meta_3D.position_to_coordinate(coords_3d[_temp], _temp))
            #         # print("Z: ", position, " - ", ides_pixel_z[i], " = ", position - ides_pixel_z[i], coords_3d[_temp] / 0.4, coordinate)
            #     else:
            #         print(coords_3d[_temp], coords_3d[_temp] / 0.4, dune_meta_3D.position_to_coordinate(coords_3d[_temp], _temp))
            #         # print("X: ", coords_3d[_temp] / 0.4, coordinate)


            value = ides_e[i]
            if abs(ides_id[i]) in track_id_map:
                cluster_index = track_id_map[abs(ides_id[i])]
            else:
                cluster_index = unknown_index
            cluster_vec_3d[cluster_index].emplace(index, value, True)


            coords_2d[0] = cos_theta*ides_z[i] + sin_theta*ides_y[i]
            coords_2d[1] = ides_x[i]

            index = dune_meta_2D[0].position_to_index(coords_2d)
            if index > dune_meta_2D[0].total_voxels():
                print("0: index ", index)
                print("   coords_2d[0]:", coords_2d[0])
                print("   coords_2d[1]:", coords_2d[1])
                print("   ides_x[i]:",    ides_x[i])
                print("   ides_y[i]:",    ides_y[i])
                print("   ides_z[i]:",    ides_z[i])
                print("   coord x: ", dune_meta_2D[0].position_to_coordinate(coords_2d[0], 0))
                print("   coord y: ", dune_meta_2D[0].position_to_coordinate(coords_2d[1], 1))
                continue

            cluster_vec_2d[0][cluster_index].emplace(index, value, True)

            coords_2d[0] = cos_theta*ides_z[i] - sin_theta*ides_y[i]
            coords_2d[1] = ides_x[i]
            index = dune_meta_2D[1].position_to_index(coords_2d)
            if index > dune_meta_2D[1].total_voxels():
                print("1: ", index)
                print("  ", coords_2d[0])
                print("  ", coords_2d[1])
                print("--", ides_y[i])
                print("--", ides_z[i])
                continue

            cluster_vec_2d[1][cluster_index].emplace(index, value, True)

            coords_2d[0] = ides_z[i]
            coords_2d[1] = ides_x[i]
            index = dune_meta_2D[2].position_to_index(coords_2d)
            if index > dune_meta_2D[2].total_voxels():
                print("2: index ", index)
                print("   coords_2d[0]:", coords_2d[0])
                print("   coords_2d[1]:", coords_2d[1])
                print("   ides_x[i]:",    ides_x[i])
                print("   ides_y[i]:",    ides_y[i])
                print("   ides_z[i]:",    ides_z[i])
                print("   coord x: ", dune_meta_2D[2].position_to_coordinate(coords_2d[0], 0))
                print("   coord y: ", dune_meta_2D[2].position_to_coordinate(coords_2d[1], 1))
                continue
            cluster_vec_2d[2][cluster_index].emplace(index, value, True)
            # if i > 10:
                # break


            #
            # if i == 1:
            #     print"(ides_x[i]: ", (ides_x[i])
            #     print"(ides_y[i]: ", (ides_y[i])
            #     print"(ides_z[i]: ", (ides_z[i])
            #     print"(coords_2d[0]: ", (coords_2d[0])
            #     print"(coords_2d[1]: ", (coords_2d[1])
            #     print"(index: ", (index)
            #     print"(dune_meta_2D[2].coordinates(index)[0]: ", (dune_meta_2D[2].coordinates(index)[0])
            #     print"(dune_meta_2D[2].coordinates(index)[1]: ", (dune_meta_2D[2].coordinates(index)[1])
            #     print"(dune_meta_2D[2].position(index)[0]: ", (dune_meta_2D[2].position(index)[0])
            #     print"(dune_meta_2D[2].position(index)[1]: ", (dune_meta_2D[2].position(index)[1])

        clusters_2d[0].emplace(cluster_vec_2d[0])
        ev_clusters2d.emplace(clusters_2d[0], dune_meta_2D[0])

        clusters_2d[1].emplace(cluster_vec_2d[1])
        ev_clusters2d.emplace(clusters_2d[1], dune_meta_2D[1])

        clusters_2d[2].emplace(cluster_vec_2d[2])
        ev_clusters2d.emplace(clusters_2d[2], dune_meta_2D[2])

        clusters_3d.emplace(cluster_vec_3d)
        ev_clusters3d.emplace(clusters_3d)


        io_manager.save_entry()

    # Save out all of the entries
    io_manager.finalize()

    return


if __name__ == '__main__':
    main()

from larcv import larcv
import numpy

def main():
    directory = "example_files/"
    # for _file in ["NC_larcv_dev.h5", "beam_larcv_dev.h5", "nueCC_larcv_dev.h5", "numuCC_larcv_dev.h5"]:
    for _file in ["anatree_beam.h5",]:
        iom = larcv.IOManager(larcv.IOManager.kREAD)
        iom.add_in_file(directory + _file)
        iom.initialize()


        for entry in range(iom.get_n_entries()):
            iom.read_entry(entry)
            neutrino  = larcv.EventParticle.to_particle(iom.get_data("particle", "neutrino"))
            particles = larcv.EventParticle.to_particle(iom.get_data("particle", "segment"))
            clusters  = larcv.EventSparseCluster2D.to_sparse_cluster(iom.get_data("cluster2d", "segment"))
            clusters3D  = larcv.EventSparseCluster3D.to_sparse_cluster(iom.get_data("cluster3d", "segment"))
            x_clusters = clusters.as_vector()[0]
            y_clusters = clusters.as_vector()[1]
            z_clusters = clusters.as_vector()[2]
            # print()
            # if (particles.size() < 25):
            #     print(entry)
            #     print(particles.size())
            # print(x_clusters.size())
            # print(y_clusters.size())
            # print(z_clusters.size())



            # for i in range(particles.size()):
            #     particle = particles.as_vector()[i]
            #     x_pix = x_clusters.as_vector()[i]
            #     y_pix = y_clusters.as_vector()[i]
            #     z_pix = z_clusters.as_vector()[i]
            #     print("ID {ID} is {pdg} by {parent}, process {proc}, E={energy} ".format(
            #         ID  = particle.track_id(), 
            #         pdg = particle.pdg_code(),
            #         parent = particle.parent_track_id(),
            #         energy = particle.energy_init(),
            #         proc   = particle.creation_process()
            #         ))

            #     print("--",x_pix.size(), numpy.sum(x_pix.values()))
            #     print("--",y_pix.size(), numpy.sum(x_pix.values()))
            #     print("--",z_pix.size(), numpy.sum(x_pix.values()))

            print(_file, neutrino.as_vector().front().pdg_code())
            print(_file, neutrino.as_vector().front().nu_current_type())

            # print(x_clusters.as_vector().back().size())
            # print(y_clusters.as_vector().back().size())
            # print(z_clusters.as_vector().back().size())
            # break


if __name__ == '__main__':
    main()
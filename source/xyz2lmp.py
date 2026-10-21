from ovito.io import import_file, export_file

node = import_file("graphene.xyz", columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])
export_file(node, "graphene.data", "lammps_data")

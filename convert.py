import re
import os
import sys
import argparse


def parse_float(value):
    """Parse a float value, handling scientific notation properly."""
    try:
        return float(value)
    except ValueError:
        # Fix incomplete scientific notation (e.g., "1.05607390e" to "1.05607390e+00")
        if value.endswith('e'):
            value += '+00'
        return float(value)


def parse_eqg(input_file):
    """Parse the EQG file into a structured format."""
    print(f"Parsing EQG file: {input_file}")

    # Storage for parsed EQG data
    eqg_data = {
        'model_name': '',
        'materials': [],
        'vertices': [],
        'faces': [],
        'bones': []
    }

    with open(input_file, 'r') as f:
        content = f.read()

    # Extract model name
    model_name_match = re.search(r'EQGMODELDEF\s+"([^"]+)"', content)
    if model_name_match:
        eqg_data['model_name'] = model_name_match.group(1)
        print(f"Model name: {eqg_data['model_name']}")

    # Extract materials
    material_section_match = re.search(
        r'NUMMATERIALS\s+(\d+)(.*?)(?=NUMVERTICES)', content, re.DOTALL)
    if material_section_match:
        num_materials = int(material_section_match.group(1))
        material_section = material_section_match.group(2)

        # Find all material blocks - improved pattern to catch the last material block
        material_blocks = re.findall(
            r'MATERIAL\s+"([^"]+)"(.*?)(?=MATERIAL\s+"|NUMVERTICES|\Z)', material_section, re.DOTALL)

        print(
            f"Found {len(material_blocks)} material blocks (expected {num_materials})")

        for material_name, block in material_blocks:
            shader_match = re.search(r'SHADERTAG\s+"([^"]+)"', block)
            shader = shader_match.group(1) if shader_match else ""

            textures = {}
            texture_matches = re.findall(
                r'PROPERTY\s+"([^"]+)"\s+\d+\s+"([^"]+)"', block)
            for prop, texture in texture_matches:
                textures[prop] = texture

            eqg_data['materials'].append({
                'name': material_name.strip(),
                'shader': shader,
                'textures': textures
            })

        print(f"Processed {len(eqg_data['materials'])} materials")

    # Verify all materials were found and processed
    if len(eqg_data['materials']) != num_materials:
        print(
            f"WARNING: Material count mismatch! Found {len(eqg_data['materials'])}, expected {num_materials}")

    # Extract vertices - IMPROVED REGEX PATTERN
    vertex_block_match = re.search(
        r'NUMVERTICES\s+(\d+)(.*?)NUMFACES', content, re.DOTALL)
    if vertex_block_match:
        num_vertices = int(vertex_block_match.group(1))
        vertex_block = vertex_block_match.group(2)

        # More robust pattern that captures the vertex index and handles the last vertex properly
        vertex_pattern = r'VERTEX\s+//\s+(\d+)([\s\S]*?)(?=VERTEX\s+//\s+\d+|NUMFACES|\Z)'
        vertex_matches = re.findall(vertex_pattern, vertex_block, re.DOTALL)

        print(
            f"Found {len(vertex_matches)} vertex blocks (expected {num_vertices})")

        for vertex_index, vertex_data in vertex_matches:
            vertex = {}
            vertex_index = int(vertex_index)

            # Extract position (XYZ)
            xyz_match = re.search(
                r'XYZ\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', vertex_data)
            if xyz_match:
                try:
                    vertex['position'] = [
                        parse_float(xyz_match.group(1)),
                        parse_float(xyz_match.group(2)),
                        parse_float(xyz_match.group(3))
                    ]
                except Exception as e:
                    print(
                        f"Warning: Could not parse XYZ for vertex {vertex_index}: {xyz_match.groups()}. Error: {e}")
                    vertex['position'] = [0.0, 0.0, 0.0]

            # Extract UV coordinates
            uv_matches = re.findall(
                r'UV\s+([-\d.e+]+)\s+([-\d.e+]+)', vertex_data)
            if uv_matches:
                vertex['uvs'] = []
                for u, v in uv_matches:
                    try:
                        vertex['uvs'].append([parse_float(u), parse_float(v)])
                    except Exception as e:
                        print(
                            f"Warning: Could not parse UV coordinates for vertex {vertex_index}: {u}, {v}. Error: {e}")
                        vertex['uvs'].append([0.0, 0.0])

            # Extract normals
            normal_match = re.search(
                r'NORMAL\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', vertex_data)
            if normal_match:
                try:
                    vertex['normal'] = [
                        parse_float(normal_match.group(1)),
                        parse_float(normal_match.group(2)),
                        parse_float(normal_match.group(3))
                    ]
                except Exception as e:
                    print(
                        f"Warning: Could not parse normal for vertex {vertex_index}: {normal_match.groups()}. Error: {e}")
                    vertex['normal'] = [0.0, 0.0, 0.0]

            # Extract color/tint
            tint_match = re.search(
                r'TINT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', vertex_data)
            if tint_match:
                vertex['tint'] = [int(tint_match.group(1)), int(tint_match.group(2)), int(
                    tint_match.group(3)), int(tint_match.group(4))]

            eqg_data['vertices'].append(vertex)

        print(f"Processed {len(eqg_data['vertices'])} vertices")

        # Verify all vertices were found and processed
        if len(eqg_data['vertices']) != num_vertices:
            print(
                f"WARNING: Vertex count mismatch! Found {len(eqg_data['vertices'])}, expected {num_vertices}")

    # Extract faces
    face_block_match = re.search(
        r'NUMFACES\s+(\d+)(.*?)(?=NUMBONES|$)', content, re.DOTALL)
    if face_block_match:
        num_faces = int(face_block_match.group(1))
        face_block = face_block_match.group(2)

        # More robust pattern for face parsing
        face_pattern = r'FACE\s+//\s+(\d+)([\s\S]*?)(?=FACE\s+//\s+\d+|NUMBONES|\Z)'
        face_matches = re.findall(face_pattern, face_block, re.DOTALL)

        print(f"Found {len(face_matches)} face blocks (expected {num_faces})")

        if len(face_matches) != num_faces:
            print("Warning: Face count mismatch. This could affect the model output.")

        for face_index, face_data in face_matches:
            face_index = int(face_index)
            face = {}

            # Extract triangle indices
            triangle_match = re.search(
                r'TRIANGLE\s+(\d+)\s+(\d+)\s+(\d+)', face_data)
            if triangle_match:
                face['indices'] = [int(triangle_match.group(1)), int(
                    triangle_match.group(3)), int(triangle_match.group(2))]
            else:
                print(
                    f"Warning: Triangle indices not found for face {face_index}")
                continue  # Skip this face if indices are missing

            # Extract material
            material_match = re.search(r'MATERIAL\s+"([^"]+)"', face_data)
            if material_match:
                face['material'] = material_match.group(1)
            else:
                # Use a default material if none is specified
                if eqg_data['materials']:
                    face['material'] = eqg_data['materials'][0]['name']
                else:
                    face['material'] = "DEFAULT"
                print(
                    f"Warning: Material not found for face {face_index}, using default: {face['material']}")

            # Extract other properties
            face['passable'] = 1 if re.search(
                r'PASSABLE\s+1', face_data) else 0
            face['transparent'] = 1 if re.search(
                r'TRANSPARENT\s+1', face_data) else 0
            face['collision'] = 1 if re.search(
                r'COLLISIONREQUIRED\s+1', face_data) else 0

            eqg_data['faces'].append(face)

        print(f"Processed {len(eqg_data['faces'])} faces")

    return eqg_data


def convert_to_s3d(eqg_data):
    """Convert the parsed EQG data to S3D format."""
    print("Converting to S3D format...")
    s3d_data = {
        'sprite_defs': [],
        'material_defs': [],
        'material_palettes': [],
        'dm_sprite_defs': [],
        'actor_def': ''
    }

    model_name = eqg_data['model_name'].upper()

    # Create sprite definitions for each material
    material_to_sprite = {}
    for i, material in enumerate(eqg_data['materials']):
        # Create a unique sprite name for this material
        sprite_name = f"{model_name}_{i}_SPRITE"
        material_to_sprite[material['name']] = sprite_name

        # Find texture for this material
        texture_file = ""
        for prop, tex in material['textures'].items():
            if prop == 'e_TextureDiffuse0':
                # Use original filename but extract basename
                texture_file = os.path.basename(tex)
                break

        # If no texture found, use a default
        if not texture_file:
            texture_file = f"{model_name}_{i}.DDS"

        # Create sprite definition
        sprite_def = f'SIMPLESPRITEDEF "{sprite_name}"\n'
        sprite_def += '\tTAGINDEX 0\n'
        sprite_def += '\tVARIATION 0\n'
        sprite_def += '\tSKIPFRAMES? NULL\n'
        sprite_def += '\tANIMATED? NULL\n'
        sprite_def += '\tSLEEP? 0\n'
        sprite_def += '\tCURRENTFRAME? NULL\n'
        sprite_def += '\tNUMFRAMES 1\n'
        sprite_def += f'\t\tFRAME "{model_name}_{i}"\n'
        sprite_def += '\t\t\tNUMFILES 1\n'
        sprite_def += f'\t\t\t\tFILE "{texture_file.upper()}"\n\n'
        s3d_data['sprite_defs'].append(sprite_def)

        # Create material definition for this sprite
        material_def = f'MATERIALDEFINITION "{model_name}_{i}_MDF"\n'
        material_def += '\tTAGINDEX 0\n'
        material_def += '\tVARIATION 0\n'
        material_def += '\tRENDERMETHOD "USERDEFINED_2"\n'
        material_def += '\tRGBPEN 178 178 178 0\n'
        material_def += '\tBRIGHTNESS 0.00000000e+00\n'
        material_def += '\tSCALEDAMBIENT 7.50000000e-01\n'  # Same as example
        material_def += '\tSIMPLESPRITEINST\n'
        material_def += f'\t\tTAG "{sprite_name}"\n'
        material_def += '\t\tSIMPLESPRITETAGINDEX 0\n'  # Added according to your example
        material_def += '\t\tHEXFIFTYFLAG 0\n'
        material_def += '\tPAIRS? 0 0.00000000e+00\n'
        material_def += '\tDOUBLESIDED 0\n\n'  # 0 instead of 1 as in your example
        s3d_data['material_defs'].append(material_def)

    # Create material palette
    material_palette = f'MATERIALPALETTE "{model_name}_MP"\n'
    material_palette += f'\tNUMMATERIALS {len(eqg_data["materials"])}\n'
    for i in range(len(eqg_data["materials"])):
        material_palette += f'\tMATERIAL "{model_name}_{i}_MDF"\n'
    material_palette += '\n'
    s3d_data['material_palettes'].append(material_palette)

    # Create DM sprite definition
    dm_sprite_def = f'DMSPRITEDEF2 "{model_name}_DMSPRITEDEF"\n'
    dm_sprite_def += '\tTAGINDEX 0\n'

    # Calculate center
    center_x = sum(v['position'][0]
                   for v in eqg_data['vertices']) / len(eqg_data['vertices'])
    center_y = sum(v['position'][1]
                   for v in eqg_data['vertices']) / len(eqg_data['vertices'])
    center_z = sum(v['position'][2]
                   for v in eqg_data['vertices']) / len(eqg_data['vertices'])

    # Set centeroffset to zero
    dm_sprite_def += f'\tCENTEROFFSET 0.00000000e+00 0.00000000e+00 0.00000000e+00\n\n'

    # Add vertices
    dm_sprite_def += f'\tNUMVERTICES {len(eqg_data["vertices"])}\n'
    for vertex in eqg_data['vertices']:
        x, y, z = vertex['position']
        dm_sprite_def += f'\t\tVXYZ {x:.8e} {y:.8e} {z:.8e}\n'

    # Add UVs
    dm_sprite_def += f'\n\tNUMUVS {len(eqg_data["vertices"])}\n'
    for vertex in eqg_data['vertices']:
        if 'uvs' in vertex and len(vertex['uvs']) > 0:
            u, v = vertex['uvs'][0]
            dm_sprite_def += f'\t\tUV {u:.8e} {v:.8e}\n'
        else:
            dm_sprite_def += '\t\tUV 0.00000000e+00 0.00000000e+00\n'

   # Add vertex normals
    dm_sprite_def += f'\n\tNUMVERTEXNORMALS {len(eqg_data["vertices"])}\n'
    for vertex in eqg_data['vertices']:
        if 'normal' in vertex:
            nx, ny, nz = vertex['normal']
            dm_sprite_def += f'\t\tNXYZ {nx:.8e} {ny:.8e} {nz:.8e}\n'
        else:
            dm_sprite_def += '\t\tNXYZ 0.00000000e+00 0.00000000e+00 0.00000000e+00\n'

    # Add vertex colors (empty section)
    dm_sprite_def += '\n\tNUMVERTEXCOLORS 0\n\n'

    # Add skin assignments (empty section)
    dm_sprite_def += '\tSKINASSIGNMENTGROUPS 0\n'

    # Add material palette reference
    dm_sprite_def += f'\tMATERIALPALETTE "{model_name}_MP"\n'
    dm_sprite_def += '\tDMTRACKINST ""\n\n'

    # Add polyhedron faces
    dm_sprite_def += '\tPOLYHEDRON\n'
    # Changed from DEFINITION to SPRITE as per example
    dm_sprite_def += '\t\tSPRITE "NEGATIVE_TWO"\n'
    dm_sprite_def += f'\tNUMFACE2S {len(eqg_data["faces"])}\n'

    for i, face in enumerate(eqg_data['faces']):
        dm_sprite_def += f'\t\tDMFACE2 //{i}\n'
        dm_sprite_def += f'\t\t\tPASSABLE {face["passable"]}\n'
        v1, v2, v3 = face['indices']
        dm_sprite_def += f'\t\t\tTRIANGLE {v1} {v2} {v3}\n'

    # Add mesh operations (empty section)
    dm_sprite_def += '\n\tNUMMESHOPS 0\n\n'

    # Create material-to-face mapping
    # Count faces per material and maintain the order they appear in the model
    material_to_faces = {}
    material_order = []

    for face in eqg_data['faces']:
        mat_name = face['material']
        if mat_name not in material_to_faces:
            material_to_faces[mat_name] = 0
            material_order.append(mat_name)
        material_to_faces[mat_name] += 1

    # Map material names to material indices in the palette
    material_name_to_index = {}
    for i, material in enumerate(eqg_data['materials']):
        material_name_to_index[material['name']] = i

    # Generate FACEMATERIALGROUPS section
    dm_sprite_def += f'\tFACEMATERIALGROUPS {len(material_order)}'
    for mat_name in material_order:
        face_count = material_to_faces[mat_name]
        mat_index = material_name_to_index.get(
            mat_name, 0)  # Default to 0 if not found
        dm_sprite_def += f' {face_count} {mat_index}'
    dm_sprite_def += '\n'

    # For VERTEXMATERIALGROUPS, we need to associate vertices with materials
    # The challenge is that vertices can be shared between faces with different materials
    # We'll use a greedy approach: assign each vertex to the first material it's used with

    # First, create a mapping of vertex to material
    vertex_to_material = {}
    for face in eqg_data['faces']:
        mat_name = face['material']
        mat_index = material_name_to_index.get(mat_name, 0)

        # Assign each vertex in this face to this material if not already assigned
        for vertex_idx in face['indices']:
            if vertex_idx not in vertex_to_material:
                vertex_to_material[vertex_idx] = mat_index

    # Count vertices per material
    material_vertex_counts = {}
    for mat_index in range(len(eqg_data['materials'])):
        material_vertex_counts[mat_index] = 0

    # Count vertices for each material
    for mat_index in vertex_to_material.values():
        material_vertex_counts[mat_index] += 1

    # Handle vertices that aren't part of any face
    unassigned_vertices = set(
        range(len(eqg_data['vertices']))) - set(vertex_to_material.keys())

    # Distribute unassigned vertices across materials proportionally
    if unassigned_vertices:
        print(
            f"Found {len(unassigned_vertices)} vertices not assigned to any face")
        if len(eqg_data['materials']) > 0:
            # Assign to first material for simplicity
            for vertex_idx in unassigned_vertices:
                vertex_to_material[vertex_idx] = 0
                material_vertex_counts[0] += 1

    # Generate VERTEXMATERIALGROUPS section
    dm_sprite_def += f'\tVERTEXMATERIALGROUPS {len(material_vertex_counts)}'
    for mat_index, count in sorted(material_vertex_counts.items()):
        if count > 0:
            dm_sprite_def += f' {count} {mat_index}'
    dm_sprite_def += '\n'

    # Calculate bounding box
    min_x = min(v['position'][0] for v in eqg_data['vertices'])
    min_y = min(v['position'][1] for v in eqg_data['vertices'])
    min_z = min(v['position'][2] for v in eqg_data['vertices'])
    max_x = max(v['position'][0] for v in eqg_data['vertices'])
    max_y = max(v['position'][1] for v in eqg_data['vertices'])
    max_z = max(v['position'][2] for v in eqg_data['vertices'])

    # Add padding to the bounding box (5% of the model's size in each dimension)
    padding_x = (max_x - min_x) * 0.05
    padding_y = (max_y - min_y) * 0.05
    padding_z = (max_z - min_z) * 0.05
    # Ensure we have at least a minimal padding of 0.001 units
    min_padding = 0.001
    padding_x = max(padding_x, min_padding)
    padding_y = max(padding_y, min_padding)
    padding_z = max(padding_z, min_padding)
    # Apply padding to bounding box
    min_x -= padding_x
    min_y -= padding_y
    min_z -= padding_z
    max_x += padding_x
    max_y += padding_y
    max_z += padding_z

    # Add bounding box
    dm_sprite_def += f'\tBOUNDINGBOXMIN {min_x:.8e} {min_y:.8e} {min_z:.8e}\n'
    dm_sprite_def += f'\tBOUNDINGBOXMAX {max_x:.8e} {max_y:.8e} {max_z:.8e}\n'

    # Calculate bounding radius (maximum distance from center to any vertex)
    max_dist = 0
    for vertex in eqg_data['vertices']:
        dx = vertex['position'][0] - center_x
        dy = vertex['position'][1] - center_y
        dz = vertex['position'][2] - center_z
        dist = (dx*dx + dy*dy + dz*dz) ** 0.5
        max_dist = max(max_dist, dist)

    # Add a 5% padding to the bounding radius
    radius_padding = max_dist * 0.05
    radius_padding = max(radius_padding, min_padding)  # Ensure minimal padding
    max_dist += radius_padding

    dm_sprite_def += f'\tBOUNDINGRADIUS {max_dist:.8e}\n\n'

    # Add flags (copy from example)
    dm_sprite_def += '\tFPSCALE 12\n'
    dm_sprite_def += '\tHEXONEFLAG 1\n'
    dm_sprite_def += '\tHEXTWOFLAG 1\n'
    dm_sprite_def += '\tHEXFOURTHOUSANDFLAG 1\n'
    dm_sprite_def += '\tHEXEIGHTTHOUSANDFLAG 0\n'
    dm_sprite_def += '\tHEXTENTHOUSANDFLAG 1\n'
    dm_sprite_def += '\tHEXTWENTYTHOUSANDFLAG 0\n\n'

    s3d_data['dm_sprite_defs'].append(dm_sprite_def)

    # Create actor definition
    actor_def = f'ACTORDEF "{model_name}_ACTORDEF"\n'
    actor_def += '\tCALLBACK "SPRITECALLBACK"\n'
    actor_def += '\tBOUNDSREF 0\n'
    actor_def += '\tCURRENTACTION? NULL\n'
    actor_def += '\tLOCATION? NULL NULL NULL NULL NULL NULL\n'
    actor_def += '\tACTIVEGEOMETRY? NULL\n'
    actor_def += '\tNUMACTIONS 1\n'
    actor_def += '\t\tACTION\n'
    actor_def += '\t\t\tUNK1 0\n'
    actor_def += '\t\t\tNUMLEVELSOFDETAIL 1\n'
    actor_def += '\t\t\t\tLEVELOFDETAIL\n'
    actor_def += f'\t\t\t\t\tSPRITE "{model_name}_DMSPRITEDEF"\n'
    actor_def += '\t\t\t\t\tSPRITEINDEX 0\n'
    actor_def += '\t\t\t\t\tMINDISTANCE 1.00000002e+30\n'
    actor_def += '\tUSEMODELCOLLIDER 0\n'  # Added from example
    actor_def += '\tUSERDATA ""\n'  # Added from example

    s3d_data['actor_def'] = actor_def

    return s3d_data


def write_s3d(s3d_data, output_file):
    """Write the converted S3D data to the output file."""
    print(f"Writing S3D file: {output_file}")

    # First, add a simple header comment
    header = "// wcemu v0.0.1\n"
    header += "// This file was created by EQG to S3D converter\n\n"

    with open(output_file, 'w') as f:
        # Write header
        f.write(header)

        # Write all sprite definitions
        for sprite_def in s3d_data['sprite_defs']:
            f.write(sprite_def)

        # Write all material definitions
        for material_def in s3d_data['material_defs']:
            f.write(material_def)

        # Write all material palettes
        for material_palette in s3d_data['material_palettes']:
            f.write(material_palette)

        # Write all DM sprite definitions
        for dm_sprite_def in s3d_data['dm_sprite_defs']:
            f.write(dm_sprite_def)

        # Write actor definition
        f.write(s3d_data['actor_def'])


def main():
    parser = argparse.ArgumentParser(
        description='Convert EQG models to S3D format.')
    parser.add_argument('input', help='Input EQG model file')
    parser.add_argument('output', help='Output S3D model file')
    parser.add_argument(
        '--it', help='Specify a custom item number for the output model name (e.g., "IT123")', default=None)

    args = parser.parse_args()

    # Run the conversion process
    eqg_data = parse_eqg(args.input)

    # If custom item number is provided, override the model name
    if args.it:
        print(f"Using custom model name: {args.it}")
        eqg_data['model_name'] = args.it

    s3d_data = convert_to_s3d(eqg_data)
    write_s3d(s3d_data, args.output)

    print("Conversion completed successfully!")


if __name__ == '__main__':
    main()

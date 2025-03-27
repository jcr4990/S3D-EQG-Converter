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


def extract_quoted_text(text):
    """Extract text from within quotes."""
    match = re.search(r'"([^"]+)"', text)
    if match:
        return match.group(1)
    return None


#############################################################
# EQG TO S3D CONVERSION FUNCTIONS
#############################################################

def parse_eqg(input_file):
    """Parse the EQG file into a structured format."""
    print(f"Parsing EQG file: {input_file}")

    with open(input_file, 'r') as f:
        content = f.read()

    # Check for multiple model definitions
    model_matches = re.findall(r'EQGMODELDEF\s+"([^"]+)"', content)

    selected_model = None
    if len(model_matches) > 1:
        print(f"Found {len(model_matches)} models in the file:")
        for i, model_name in enumerate(model_matches):
            print(f"  {i+1}. {model_name}")

        while True:
            try:
                selection = int(
                    input("Enter the number of the model you want to parse and convert: "))
                if 1 <= selection <= len(model_matches):
                    selected_model = model_matches[selection-1]
                    print(f"Selected model: {selected_model}")
                    break
                else:
                    print(
                        f"Please enter a number between 1 and {len(model_matches)}")
            except ValueError:
                print("Please enter a valid number")
    else:
        if model_matches:
            selected_model = model_matches[0]
            print(f"Found single model: {selected_model}")
        else:
            print("No models found in the file!")
            return None

    # If we found multiple models, split the content to only process the selected one
    if len(model_matches) > 1:
        # Find the start position of the selected model definition
        model_start_pattern = f'EQGMODELDEF\\s+"{selected_model}"'
        start_match = re.search(model_start_pattern, content)
        if not start_match:
            print(
                f"ERROR: Could not locate the selected model {selected_model} in the file!")
            return None

        start_pos = start_match.start()

        # Find the start position of the next model definition (if any)
        next_model_pattern = r'EQGMODELDEF\s+"[^"]+"'
        next_matches = list(re.finditer(
            next_model_pattern, content[start_pos+1:]))

        if next_matches:
            # Extract content from start of selected model to start of next model
            end_pos = start_pos + 1 + next_matches[0].start()
            model_content = content[start_pos:end_pos]
        else:
            # Extract content from start of selected model to end of file
            model_content = content[start_pos:]
    else:
        # Process the entire file
        model_content = content

    # Storage for parsed EQG data
    eqg_data = {
        'model_name': selected_model,
        'materials': [],
        'vertices': [],
        'faces': [],
        'bones': []
    }

    print(f"Model name: {eqg_data['model_name']}")

    # Extract materials
    material_section_match = re.search(
        r'NUMMATERIALS\s+(\d+)(.*?)(?=NUMVERTICES)', model_content, re.DOTALL)
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
        r'NUMVERTICES\s+(\d+)(.*?)NUMFACES', model_content, re.DOTALL)
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
        r'NUMFACES\s+(\d+)(.*?)(?=NUMBONES|$)', model_content, re.DOTALL)
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

    # Extract bones
    bone_block_match = re.search(
        r'NUMBONES\s+(\d+)(.*?)(?=\Z)', model_content, re.DOTALL)
    if bone_block_match:
        num_bones = int(bone_block_match.group(1))
        bone_block = bone_block_match.group(2)

        # Pattern to match bone definitions
        bone_pattern = r'BONE\s+"([^"]+)"\s*//\s*(\d+)([\s\S]*?)(?=BONE\s+"|$)'
        bone_matches = re.findall(bone_pattern, bone_block, re.DOTALL)

        print(f"Found {len(bone_matches)} bone blocks (expected {num_bones})")

        for bone_name, bone_index, bone_data in bone_matches:
            bone = {
                'name': bone_name,
                'index': int(bone_index)
            }

            # Extract next bone
            next_match = re.search(r'NEXT\s+(-?\d+)', bone_data)
            if next_match:
                bone['next'] = int(next_match.group(1))

            # Extract children
            children_match = re.search(r'CHILDREN\s+(\d+)', bone_data)
            if children_match:
                bone['children'] = int(children_match.group(1))

            # Extract child index
            child_index_match = re.search(r'CHILDINDEX\s+(-?\d+)', bone_data)
            if child_index_match:
                bone['child_index'] = int(child_index_match.group(1))

            # Extract pivot
            pivot_match = re.search(
                r'PIVOT\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', bone_data)
            if pivot_match:
                try:
                    bone['pivot'] = [
                        parse_float(pivot_match.group(1)),
                        parse_float(pivot_match.group(2)),
                        parse_float(pivot_match.group(3))
                    ]
                except Exception as e:
                    print(
                        f"Warning: Could not parse pivot for bone {bone_name}: {e}")
                    bone['pivot'] = [0.0, 0.0, 0.0]

            # Extract quaternion
            quat_match = re.search(
                r'QUATERNION\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', bone_data)
            if quat_match:
                try:
                    bone['quaternion'] = [
                        parse_float(quat_match.group(1)),
                        parse_float(quat_match.group(2)),
                        parse_float(quat_match.group(3)),
                        parse_float(quat_match.group(4))
                    ]
                except Exception as e:
                    print(
                        f"Warning: Could not parse quaternion for bone {bone_name}: {e}")
                    bone['quaternion'] = [0.0, 0.0, 0.0, 1.0]

            # Extract scale
            scale_match = re.search(
                r'SCALE\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', bone_data)
            if scale_match:
                try:
                    bone['scale'] = [
                        parse_float(scale_match.group(1)),
                        parse_float(scale_match.group(2)),
                        parse_float(scale_match.group(3))
                    ]
                except Exception as e:
                    print(
                        f"Warning: Could not parse scale for bone {bone_name}: {e}")
                    bone['scale'] = [1.0, 1.0, 1.0]

            eqg_data['bones'].append(bone)

        print(f"Processed {len(eqg_data['bones'])} bones")

        # Verify all bones were found and processed
        if len(eqg_data['bones']) != num_bones:
            print(
                f"WARNING: Bone count mismatch! Found {len(eqg_data['bones'])}, expected {num_bones}")

    return eqg_data


def convert_to_s3d(eqg_data):
    """Convert the parsed EQG data to S3D format."""
    if not eqg_data:
        print("No EQG data provided for conversion. Aborting.")
        return None

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
    # Determine if we need to flip the UVs
    need_flip = False

    # Auto-detect by analyzing the UVs in the model
    if len(eqg_data["vertices"]) > 10:
        # Check a sample of vertices with UVs
        sample_vertices = []
        for v in eqg_data["vertices"]:
            if 'uvs' in v and len(v['uvs']) > 0:
                sample_vertices.append(v)
                if len(sample_vertices) >= 20:  # Sample size of 20 for better accuracy
                    break

        if sample_vertices:
            # Look at distribution of UV values
            # In many asset creation tools, UV Y values are flipped (1 at bottom rather than top)
            v_values = [v['uvs'][0][1] for v in sample_vertices]
            v_in_range = sum(1 for v in v_values if 0 <= v <= 1)

            # If most are already in 0-1 range, they likely need flipping
            if v_in_range >= len(v_values) * 0.7:  # 70% threshold
                need_flip = True

            print(
                f"UV analysis: {v_in_range}/{len(v_values)} Y values in 0-1 range. {'Flipping' if need_flip else 'Not flipping'} UVs.")

    dm_sprite_def += f'\n\tNUMUVS {len(eqg_data["vertices"])}\n'
    for vertex in eqg_data['vertices']:
        if 'uvs' in vertex and len(vertex['uvs']) > 0:
            u, v = vertex['uvs'][0]
            # Flip the V coordinate (1.0 - v) if needed
            if need_flip:
                v = 1.0 - v
            dm_sprite_def += f'\t\tUV {u:.8e} {v:.8e}\n'
        else:
            # Default UV, potentially y-flipped
            v_default = 1.0 if need_flip else 0.0
            dm_sprite_def += f'\t\tUV 0.00000000e+00 {v_default:.8e}\n'

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
    material_face_groups = []
    for mat_name in material_order:
        face_count = material_to_faces[mat_name]
        mat_index = material_name_to_index.get(
            mat_name, 0)  # Default to 0 if not found
        material_face_groups.append((face_count, mat_index))

    dm_sprite_def += f'\tFACEMATERIALGROUPS {len(material_face_groups)}'
    for face_count, mat_index in material_face_groups:
        dm_sprite_def += f' {face_count} {mat_index}'
    dm_sprite_def += '\n'

    # For VERTEXMATERIALGROUPS, we need to associate vertices with materials
    # The challenge is that vertices can be shared between faces with different materials
    # We'll use a greedy approach: assign each vertex to the first material it's used with

    # First, create a mapping of vertex to material
    vertex_to_material = {}
    for face in eqg_data['faces']:
        mat_name = face['material']
        if mat_name in material_name_to_index:
            mat_index = material_name_to_index[mat_name]

            # Assign each vertex in this face to this material if not already assigned
            for vertex_idx in face['indices']:
                if vertex_idx not in vertex_to_material and vertex_idx < len(eqg_data['vertices']):
                    vertex_to_material[vertex_idx] = mat_index

    # Count vertices per material
    material_vertex_counts = {}
    for mat_index in range(len(eqg_data['materials'])):
        material_vertex_counts[mat_index] = 0

    # Count vertices for each material
    for vertex_idx, mat_index in vertex_to_material.items():
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
    material_vertex_groups = []
    for mat_index, count in material_vertex_counts.items():
        if count > 0:
            material_vertex_groups.append((count, mat_index))

    dm_sprite_def += f'\tVERTEXMATERIALGROUPS {len(material_vertex_groups)}'
    for count, mat_index in material_vertex_groups:
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
    header = "// This Quail .wce file was created by EQG to S3D converter\n\n"

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


#############################################################
# S3D TO EQG CONVERSION FUNCTIONS
#############################################################

def parse_s3d(input_file):
    """Parse the S3D file into a structured format."""
    print(f"Parsing S3D file: {input_file}")

    try:
        with open(input_file, 'r') as f:
            content = f.read()

        # Find actor definition to get model name
        actor_def_match = re.search(r'ACTORDEF\s+"([^"]+)"', content)
        if not actor_def_match:
            print("No ACTORDEF found in the file! Cannot determine model name.")
            return None

        actor_def_name = actor_def_match.group(1)
        # Extract model name from ACTORDEF name (usually IT#_ACTORDEF)
        model_name = actor_def_name.replace("_ACTORDEF", "")
        print(f"Found model name: {model_name}")

        # Find the DMSPRITEDEF2 block
        dm_sprite_match = re.search(
            r'DMSPRITEDEF2\s+"([^"]+)"(.*?)(?=ACTORDEF|\Z)', content, re.DOTALL)
        if not dm_sprite_match:
            print("No DMSPRITEDEF2 found in the file! Cannot extract model data.")
            return None

        dm_sprite_name = dm_sprite_match.group(1)
        dm_sprite_content = dm_sprite_match.group(2)
        print(f"Found DMSPRITEDEF2: {dm_sprite_name}")

        # Find material palette
        material_palette_match = re.search(
            r'MATERIALPALETTE\s+"([^"]+)"(.*?)(?=DMSPRITEDEF2|\Z)', content, re.DOTALL)
        if not material_palette_match:
            print("No MATERIALPALETTE found in the file! Cannot extract materials.")
            return None

        material_palette_name = material_palette_match.group(1)
        material_palette_content = material_palette_match.group(2)
        print(f"Found MATERIALPALETTE: {material_palette_name}")

        # Extract materials
        material_defs = []
        material_matches = re.findall(
            r'MATERIAL\s+"([^"]+)"', material_palette_content)
        print(f"Found {len(material_matches)} materials in palette")

        # For each material in the palette, find its definition
        for material_name in material_matches:
            material_def_match = re.search(
                f'MATERIALDEFINITION\\s+"{material_name}"(.*?)(?=MATERIALDEFINITION|SIMPLESPRITEDEF|MATERIALPALETTE|DMSPRITEDEF2|\\Z)', content, re.DOTALL)
            if material_def_match:
                material_def_content = material_def_match.group(1)

                # Find the sprite tag
                sprite_tag_match = re.search(
                    r'TAG\s+"([^"]+)"', material_def_content)
                sprite_tag = sprite_tag_match.group(
                    1) if sprite_tag_match else ""

                # Find the sprite definition
                sprite_def_match = re.search(
                    f'SIMPLESPRITEDEF\\s+"{sprite_tag}"(.*?)(?=SIMPLESPRITEDEF|MATERIALDEFINITION|\\Z)', content, re.DOTALL)
                if sprite_def_match:
                    sprite_def_content = sprite_def_match.group(1)

                    # Find the texture file
                    file_match = re.search(
                        r'FILE\s+"([^"]+)"', sprite_def_content)
                    texture_file = file_match.group(1) if file_match else ""

                    material_defs.append({
                        'name': material_name,
                        'sprite_tag': sprite_tag,
                        'texture_file': texture_file
                    })
                    print(
                        f"  Processed material: {material_name}, texture: {texture_file}")

        # Create data structure for EQG model
        s3d_data = {
            'model_name': model_name,
            'materials': material_defs,
            'vertices': [],
            'faces': [],
            'bones': []  # S3D may not have bone data
        }

        # Extract vertices
        vertices_match = re.search(
            r'NUMVERTICES\s+(\d+)(.*?)NUMUVS', dm_sprite_content, re.DOTALL)
        if vertices_match:
            num_vertices = int(vertices_match.group(1))
            vertices_content = vertices_match.group(2)
            print(f"Found NUMVERTICES: {num_vertices}")

            vertex_pattern = r'VXYZ\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)'
            vertex_matches = re.findall(vertex_pattern, vertices_content)
            print(f"Extracted {len(vertex_matches)} vertex positions")

            for i, (x, y, z) in enumerate(vertex_matches):
                s3d_data['vertices'].append({
                    'position': [parse_float(x), parse_float(y), parse_float(z)],
                    'uvs': [],
                    'normal': []
                })

        # Extract UVs
        uvs_match = re.search(
            r'NUMUVS\s+(\d+)(.*?)NUMVERTEXNORMALS', dm_sprite_content, re.DOTALL)
        if uvs_match:
            num_uvs = int(uvs_match.group(1))
            uvs_content = uvs_match.group(2)
            print(f"Found NUMUVS: {num_uvs}")

            uv_pattern = r'UV\s+([-\d.e+]+)\s+([-\d.e+]+)'
            uv_matches = re.findall(uv_pattern, uvs_content)
            print(f"Extracted {len(uv_matches)} UV coordinates")

            for i, (u, v) in enumerate(uv_matches):
                if i < len(s3d_data['vertices']):
                    # Preserve the original sign of the V coordinate
                    u_val = parse_float(u)
                    v_val = parse_float(v)
                    s3d_data['vertices'][i]['uvs'].append([u_val, v_val])

        # Extract normals
        normals_match = re.search(
            r'NUMVERTEXNORMALS\s+(\d+)(.*?)(?=NUMVERTEXCOLORS|SKINASSIGNMENTGROUPS)', dm_sprite_content, re.DOTALL)
        if normals_match:
            num_normals = int(normals_match.group(1))
            normals_content = normals_match.group(2)
            print(f"Found NUMVERTEXNORMALS: {num_normals}")

            normal_pattern = r'NXYZ\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)'
            normal_matches = re.findall(normal_pattern, normals_content)
            print(f"Extracted {len(normal_matches)} vertex normals")

            for i, (nx, ny, nz) in enumerate(normal_matches):
                if i < len(s3d_data['vertices']):
                    s3d_data['vertices'][i]['normal'] = [
                        parse_float(nx), parse_float(ny), parse_float(nz)]

        # Extract faces
        faces_match = re.search(
            r'NUMFACE2S\s+(\d+)(.*?)(?=NUMMESHOPS|\Z)', dm_sprite_content, re.DOTALL)
        if faces_match:
            num_faces = int(faces_match.group(1))
            faces_content = faces_match.group(2)
            print(f"Found NUMFACE2S: {num_faces}")

            # Extract face material groups to assign materials to faces
            material_assignments = []

            material_groups_match = re.search(
                r'FACEMATERIALGROUPS\s+(\d+)(.*?)(?=VERTEXMATERIALGROUPS|\Z)', dm_sprite_content, re.DOTALL)
            if material_groups_match:
                num_groups = int(material_groups_match.group(1))
                groups_content = material_groups_match.group(2).strip()
                print(f"Found FACEMATERIALGROUPS with {num_groups} groups")
                print(f"Raw material groups content: {groups_content}")

                # Extract only the numeric values
                numbers = re.findall(r'\d+', groups_content)
                print(f"Extracted numbers: {numbers}")

                if len(numbers) >= 2:  # Ensure we have at least one pair of numbers
                    try:
                        i = 0
                        while i < len(numbers) - 1:
                            face_count = int(numbers[i])
                            material_index = int(numbers[i+1])

                            for _ in range(face_count):
                                material_assignments.append(material_index)

                            i += 2
                    except Exception as e:
                        print(f"Error parsing FACEMATERIALGROUPS: {str(e)}")
                        # Create a fallback assignment - assign all faces to material 0
                        material_assignments = [0] * num_faces
                        print(
                            "Using fallback material assignments (all faces assigned to material 0)")
                else:
                    # Not enough numbers, create a fallback
                    material_assignments = [0] * num_faces
                    print(
                        "Not enough numbers in FACEMATERIALGROUPS, using fallback assignments")
            else:
                # No material groups found, create a fallback
                material_assignments = [0] * num_faces
                print("No FACEMATERIALGROUPS found, using fallback assignments")

            print(
                f"Generated {len(material_assignments)} material assignments for faces")

            # Process each face
            # More robust pattern that doesn't assume TRIANGLE directly follows PASSABLE
            face_pattern = r'DMFACE2\s+//(\d+).*?PASSABLE\s+(\d+).*?TRIANGLE\s+(\d+)\s+(\d+)\s+(\d+)'
            face_matches = re.findall(face_pattern, faces_content, re.DOTALL)
            print(f"Extracted {len(face_matches)} faces")

            for i, (face_idx, passable, v1, v2, v3) in enumerate(face_matches):
                face_idx = int(face_idx)
                material_index = material_assignments[face_idx] if face_idx < len(
                    material_assignments) else 0

                material_name = material_defs[material_index]['name'] if material_index < len(
                    material_defs) else "DEFAULT"

                # Convert the triangle indices - may need to adjust ordering from S3D to EQG
                s3d_data['faces'].append({
                    # Adjusted from S3D to EQG winding order
                    'indices': [int(v1), int(v3), int(v2)],
                    'material': material_name.replace("_MDF", ""),
                    'passable': int(passable),
                    'transparent': 0,  # Default values as these may not be directly available in S3D
                    'collision': 0
                })

                # Print the first and last 5 faces for debugging
                if i < 5 or i >= len(face_matches) - 5:
                    print(
                        f"  Face {face_idx}: material={material_name}, vertices={v1},{v2},{v3}, passable={passable}")

        print(
            f"Completed parsing S3D model with {len(s3d_data['vertices'])} vertices, {len(s3d_data['faces'])} faces, {len(s3d_data['materials'])} materials")
        return s3d_data

    except Exception as e:
        print(f"Error parsing S3D file: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def convert_to_eqg(s3d_data):
    """Convert the parsed S3D data to EQG format."""
    if not s3d_data:
        print("No S3D data provided for conversion. Aborting.")
        return None

    print("Converting to EQG format...")

    eqg_output = []

    # Start with EQGMODELDEF header
    eqg_output.append(f'EQGMODELDEF "{s3d_data["model_name"]}"')
    eqg_output.append('\tVERSION 3')

    # Add materials
    eqg_output.append(f'\tNUMMATERIALS {len(s3d_data["materials"])}')
    for i, material in enumerate(s3d_data["materials"]):
        material_name = material["name"].replace("_MDF", "")
        eqg_output.append(f'\t\tMATERIAL "{material_name}"')
        eqg_output.append('\t\t\tSHADERTAG "Opaque_MaxCB1.fx"')
        eqg_output.append('\t\t\tHEXONEFLAG 0')
        eqg_output.append('\t\t\tNUMPROPERTIES 1')

        # Add texture property
        texture_file = material["texture_file"].lower()
        eqg_output.append(
            f'\t\t\t\tPROPERTY "e_TextureDiffuse0" 2 "{texture_file}"')

        # Add animation properties
        eqg_output.append('\t\t\tANIMSLEEP 0')
        eqg_output.append('\t\t\tNUMANIMTEXTURES 0')

    # Add vertices
    eqg_output.append(f'\tNUMVERTICES {len(s3d_data["vertices"])}')
    for i, vertex in enumerate(s3d_data["vertices"]):
        eqg_output.append(f'\t\tVERTEX // {i}')

        # Position
        x, y, z = vertex["position"]
        eqg_output.append(f'\t\t\tXYZ {x:.8e} {y:.8e} {z:.8e}')

        # UVs
        if vertex["uvs"]:
            u, v = vertex["uvs"][0]
            eqg_output.append(f'\t\t\tUV {u:.8e} {v:.8e}')
            eqg_output.append('\t\t\tUV2 0.00000000e+00 0.00000000e+00')
        else:
            eqg_output.append('\t\t\tUV 0.00000000e+00 0.00000000e+00')
            eqg_output.append('\t\t\tUV2 0.00000000e+00 0.00000000e+00')

        # Normal
        if vertex["normal"]:
            nx, ny, nz = vertex["normal"]
            eqg_output.append(f'\t\t\tNORMAL {nx:.8e} {ny:.8e} {nz:.8e}')
        else:
            eqg_output.append(
                '\t\t\tNORMAL 0.00000000e+00 1.00000000e+00 0.00000000e+00')

        # Add vertex color (default)
        eqg_output.append('\t\t\tTINT 120 120 120 255')
        eqg_output.append('\t\t\tNUMWEIGHTS 0')

    # Add faces
    eqg_output.append(f'\tNUMFACES {len(s3d_data["faces"])}')
    for i, face in enumerate(s3d_data["faces"]):
        eqg_output.append(f'\t\tFACE // {i}')

        # Add triangle indices
        v1, v2, v3 = face["indices"]
        eqg_output.append(f'\t\t\tTRIANGLE {v1} {v2} {v3}')

        # Add material and flags
        eqg_output.append(f'\t\t\tMATERIAL "{face["material"]}"')
        eqg_output.append(f'\t\t\tPASSABLE {face["passable"]}')
        eqg_output.append(f'\t\t\tTRANSPARENT {face["transparent"]}')
        eqg_output.append(f'\t\t\tCOLLISIONREQUIRED {face["collision"]}')
        eqg_output.append('\t\t\tCULLED 0')
        eqg_output.append('\t\t\tDEGENERATE 0')

    # Add bones section (even if empty)
    eqg_output.append('\tNUMBONES 0')

    return '\n'.join(eqg_output)


def write_eqg(eqg_content, output_file):
    """Write the converted EQG data to the output file."""
    print(f"Writing EQG file: {output_file}")

    with open(output_file, 'w') as f:
        f.write(eqg_content)


def main():
    parser = argparse.ArgumentParser(
        description='Convert between EQG and S3D model formats.')
    parser.add_argument('input', help='Input model file')
    parser.add_argument('output', help='Output model file')
    parser.add_argument(
        '--it', help='Specify a custom item number/name for the output model', default=None)

    args = parser.parse_args()
    direction = ""

    with open(args.input, 'r') as f:
        # Read first 500 chars to determine file type
        content = f.read(500)
        if 'EQGMODELDEF' in content:
            direction = 'to_s3d'
        elif 'SIMPLESPRITEDEF' in content or 'MATERIALDEFINITION' in content:
            direction = 'to_eqg'
        else:
            print("Could not determine file type.")
            return 1

    # Run the appropriate conversion process
    if direction == 'to_s3d':
        print("Converting EQG to S3D...")
        # Parse EQG
        eqg_data = parse_eqg(args.input)
        if not eqg_data:
            print("ERROR: Failed to parse EQG data. Conversion aborted.")
            return 1

        # Apply custom name if provided
        if args.it:
            print(f"Using custom model name: {args.it}")
            eqg_data['model_name'] = args.it

        # Convert to S3D
        s3d_data = convert_to_s3d(eqg_data)
        if not s3d_data:
            print("ERROR: Failed to convert EQG data to S3D format. Conversion aborted.")
            return 1

        # Write output
        write_s3d(s3d_data, args.output)

    elif direction == 'to_eqg':
        print("Converting S3D to EQG...")
        # Parse S3D
        s3d_data = parse_s3d(args.input)
        if not s3d_data:
            print("ERROR: Failed to parse S3D data. Conversion aborted.")
            return 1

        # Apply custom name if provided
        if args.it:
            print(f"Using custom model name: {args.it}")
            s3d_data['model_name'] = args.it

        # Convert to EQG
        eqg_content = convert_to_eqg(s3d_data)
        if not eqg_content:
            print("ERROR: Failed to convert S3D data to EQG format. Conversion aborted.")
            return 1

        # Write output
        write_eqg(eqg_content, args.output)

    print("Conversion completed successfully!")
    return 0


if __name__ == '__main__':
    main()

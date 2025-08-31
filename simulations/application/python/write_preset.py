import json

def parse_freesurfer_lut(lut_path):
    annotations = []
    indexed_colors = []

    with open(lut_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                label_id = int(parts[0])
                label_name = parts[1]
                r, g, b = map(int, parts[2:5])

                annotations.append(str(label_id))
                annotations.append(label_name)

                indexed_colors.extend([
                    r / 255.0,
                    g / 255.0,
                    b / 255.0
                ])
            except ValueError:
                continue

    return annotations, indexed_colors

def build_paraview_preset(name, annotations, indexed_colors):
    preset = {
        "Name": name,
        "ColorSpace": "RGB",
        "IndexedColors": indexed_colors,
        "Annotations": annotations,
        "UseDiscreteColors": True,
        "NanColor": [0.5, 0.5, 0.5]
    }
    return [preset]  # wrap in list!

# === Modify this path ===
lut_path = "/usr/local/freesurfer/7.4.1/FreeSurferColorLUT.txt"

annotations, indexed_colors = parse_freesurfer_lut(lut_path)
preset_list = build_paraview_preset("FreeSurferLUT", annotations, indexed_colors)

with open("freesurfer_preset.json", "w") as f:
    json.dump(preset_list, f, indent=2)

print("ParaView-compatible preset saved to 'freesurfer_preset.json'")


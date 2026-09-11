from __future__ import annotations
import numpy as np
import theia.units as u
import theia.material
from material_ice import LayerMediumModel
from theia.surface import BorderSurface
from theia.surface import AbsorbingSurface
from theia.material import Material
from theia.scene import MeshStore
import trimesh
from theia.scene import _createMeshFromTrimesh
from theia.material import MaterialStore
from theia.scene import Transform
from theia.scene import Scene
from theia.volume import Attenuating
from theia.material import MaterialFlags


class LayerStack:

    def __init__(self, filename):
        
        self.filename = filename
        
        self.layer_thickness = 10 * u.m
        self.width=1000*u.m
        self.length=1000*u.m
        
        self.g = 0.9
        self.n = 1.15
        self.ng = 1.35

        self.media = []
        self.materials = []
        self.instances = []

        self._load_file()
        self._create_media()
        self._create_materials()
        self._create_geometry()
        self.scene = self.build_scene()

        
    def _load_file(self):
    
        data = np.loadtxt(self.filename)
    
        self.depth = data[:, 0]
        self.z = 1948.07 - self.depth
        
        self.be400 = data[:, 1]
        self.adust400 = data[:, 2]
        self.delta_tau = data[:, 3]


    def _create_media(self):

        self.media = []
    
        for i in range(len(self.depth)):
    
            model = LayerMediumModel(
                be400=self.be400[i],
                adust400=self.adust400[i],
                delta_tau=self.delta_tau[i],
                name=f"Layer_{i}",
            )
    
            medium = model.createMedium(
                physicModel=Attenuating(),
                numSamples=8192,
                #numSamples=16384,
            )
    
            self.media.append(medium)
        
    def _create_materials(self):
    
        self.materials = []
    
        for i in range(len(self.media)):
    
            if i == 0:
                outside = None # draußen ist Vakuum
            else:
                outside = self.media[i - 1]
    
            material = Material(
                name=f"Boundary_{i}",
                inside=self.media[i],
                outside=self.media[i - 1],
                physicModel=BorderSurface(),
                flags=MaterialFlags.VOLUME_BORDER,
            )
            self.materials.append(material)
    
        # DOM
        dom_layer = self.get_layer(-383.40)
    
        dom_material = Material(
            name="DOM",
            inside=self.media[dom_layer],
            outside=self.media[dom_layer],
            physicModel=AbsorbingSurface(),
            flags=(
                MaterialFlags.DETECTOR
                | MaterialFlags.BLACK_BODY
                | MaterialFlags.SKIP_MEDIA_MISMATCH_TEST
            ),
        )
    
        self.materials.append(dom_material)
    
        self.material_store = MaterialStore(
            material=self.materials
        )
                
    def _create_geometry(self):
    
        half_width = 0.5 * self.width
        half_length = 0.5 * self.length

        #Hier sind die Ecken der horizonateln Fläche bei z=0
        vertices = np.array([
            [-half_width, -half_length, 0.0],
            [ half_width, -half_length, 0.0],
            [ half_width,  half_length, 0.0],
            [-half_width,  half_length, 0.0],
        ])
        #triangle faces (2 Dreiecke)
        faces = np.array([
            [0, 1, 2],
            [0, 2, 3],
        ])
        
        plane = trimesh.Trimesh(
            vertices=vertices,
            faces=faces,
        )
        
        plane_mesh = _createMeshFromTrimesh(plane)
        
        
        sphere = trimesh.creation.icosphere(
            subdivisions=3,
            radius=0.5,
        )
        
        sphere_mesh = _createMeshFromTrimesh(sphere)
        
    
        self.mesh_store = MeshStore({
            "boundary": plane_mesh,
            "dom": sphere_mesh,
        })

        #Jetzt dieselbe plane mit verschiednenen Translationen zu instances hinzufügen 
        self.instances = []
        
        for i in range(1, len(self.media)):
        
            boundary_z = 0.5 * (
                self.z[i - 1] + self.z[i]
            )
        
            instance = self.mesh_store.createInstance(
                key="boundary",
                material=f"Boundary_{i}",
                transform=Transform.Translation(
                    0.0, 0.0, boundary_z
                ),
            )
            
            self.instances.append(instance)
        
        
    
        dom = self.mesh_store.createInstance(
            key="dom",
            material="DOM",
            transform=Transform.Translation(
                31.25,
                -72.93,
                -383.40,
            ),
            detectorId=1,
        )
        self.instances.append(dom)
        
    def build_scene(self):

        self.scene = Scene(
            instances=self.instances,
            materials=self.material_store,
        )
        
        return self.scene

    def get_layer(self, z):

        idx = np.argmin(np.abs(self.z - z))
        return idx

    def get_medium_index(self, z):

        idx = self.get_layer(z)
        return self.material_store.media[self.media[idx].name]
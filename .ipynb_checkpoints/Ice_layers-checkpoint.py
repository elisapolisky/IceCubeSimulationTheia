from __future__ import annotations
import numpy as np
import theia.units as u
import theia.material
from LayerMediumModel import LayerMediumModel
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
        self.mu_s = data[:, 1]
        self.mu_a = data[:, 2]
        #self.be400 = data[:, 1]
        #self.adust400 = data[:, 2]
        #self.delta_tau = data[:, 3]

    def _create_media(self):
        self.media = []
        for i in range(len(self.depth)):
            model = LayerMediumModel(
                mu_a=self.mu_a[i],
                mu_s=self.mu_s[i],
                g=self.g,
                n=self.n,
                ng=self.ng,
                name=f"Layer_{i}",
            )
            medium = model.createMedium(physicModel = Attenuating())
            self.media.append(medium)

    #def _create_media(self):

        #self.media = []
    
        #for i in range(len(self.depth)):
    
            #model = LayerMediumModel(
                #be400=self.be400[i],
                #adust400=self.adust400[i],
                #delta_tau=self.delta_tau[i],
                #name=f"Layer_{i}",
            #)
    
            #medium = model.createMedium(
             #   physicModel=Attenuating()
            #)
    
            #self.media.append(medium)
        
    def _create_materials(self):

        self.materials = []
    
        for i, medium in enumerate(self.media):
    
            material = Material(
                name=f"Layer_{i}",
                inside=medium,
                outside=medium,
                physicModel= BorderSurface(),
                #flags=...
            )

            self.materials.append(material)
            
        dom_layer = self.get_layer(-383.40)    
        dom_material = Material(
            name="DOM",
            inside=self.media[dom_layer],
            outside=self.media[dom_layer],
            physicModel=AbsorbingSurface(),
            flags= (MaterialFlags.DETECTOR
                    | MaterialFlags.BLACK_BODY
                    | MaterialFlags.SKIP_MEDIA_MISMATCH_TEST
            )
        )
        self.materials.append(dom_material)
        
        self.material_store = MaterialStore(
            material=self.materials
        )
            
    def _create_geometry(self):

        box = trimesh.creation.box(
            extents=(
                self.width,
                self.length,
                self.layer_thickness,
            )
        )
        sphere = trimesh.creation.icosphere(
            subdivisions=2,
            radius=0.5,
        )
            
        box_mesh = _createMeshFromTrimesh(box)
        sphere_mesh = _createMeshFromTrimesh(sphere)
    
        self.mesh_store = MeshStore({
            "layer": box_mesh,
            "dom": sphere_mesh,
        })
    
        self.instances = []
    
        for i in range(len(self.depth)):
    
            transform = Transform.Translation(
                0.0,
                0.0,
                self.z[i]
            )
    
            instance = self.mesh_store.createInstance(
                key="layer",
                material=f"Layer_{i}",
                transform=transform,
            )
        
        
            self.instances.append(instance)
        
        layer = self.get_layer(-383.40)
        dom = self.mesh_store.createInstance(
            key="dom",
            material= "DOM",      
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
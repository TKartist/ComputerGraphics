// definition of a plane
var plane_vertices = [
    -0.5,  0.0,  -0.5,
    -0.5, 0.0,  0.5,
     0.5, 0.0,  0.5,
    -0.5,  0.0,  -0.5,
     0.5, 0.0,  0.5,
     0.5,  0.0,  -0.5,
];
var plane_normals = [
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
];

var plane_colors = [
    0.56, 0.45, 0.4,
    0.56, 0.45, 0.4,
    0.56, 0.45, 0.4,
    0.56, 0.45, 0.4,
    0.56, 0.45, 0.4,
    0.56, 0.45, 0.4,
];

//---------------------------
// definition of the cube
//---------------------------
var cube_vertices = [
    -0.5,  0.5,  -0.5, // front face
    -0.5, -0.5,  -0.5,
     0.5, -0.5,  -0.5,
    -0.5,  0.5,  -0.5,
     0.5, -0.5,  -0.5,
     0.5,  0.5,  -0.5,
     0.5, -0.5,  0.5,// back face
     -0.5, -0.5,  0.5,
    -0.5,  0.5,  0.5,
     0.5,  0.5,  0.5,
     0.5, -0.5,  0.5,
    -0.5,  0.5,  0.5,
     0.5, -0.5, -0.5,// right face
     0.5, -0.5,  0.5,
     0.5,  0.5,  0.5,
     0.5, 0.5, -0.5,
     0.5, -0.5, -0.5,
     0.5, 0.5, 0.5,
     -0.5,  0.5,  0.5,// left face
     -0.5, -0.5,  0.5,
     -0.5, -0.5, -0.5,
     -0.5, 0.5, 0.5,
     -0.5, -0.5, -0.5,
     -0.5, 0.5, -0.5,
     -0.5, 0.5, -0.5, //top face
     0.5, 0.5, 0.5,
     -0.5, 0.5, 0.5,
     -0.5, 0.5, -0.5,
     0.5, 0.5, -0.5,
     0.5, 0.5, 0.5,
     -0.5, -0.5, -0.5, //bottom face
     -0.5, -0.5, 0.5,
     0.5, -0.5, 0.5,
     -0.5, -0.5, -0.5,
     0.5, -0.5, 0.5,
     0.5, -0.5, -0.5,
];

var cube_colors = [
    0.9,  0.24,  0.24, // front face
    0.9,  0.24,  0.24,
    0.9,  0.24,  0.24,
    0.9,  0.24,  0.24,
    0.9,  0.24,  0.24,
    0.9,  0.24,  0.24,
    0.6,  0.15,  0.71, // back face
    0.6,  0.15,  0.71,
    0.6,  0.15,  0.71,
    0.6,  0.15,  0.71,
    0.6,  0.15,  0.71,
    0.6,  0.15,  0.71,
    0.0, 0.61, 0.87, // right face
    0.0, 0.61, 0.87,
    0.0, 0.61, 0.87,
    0.0, 0.61, 0.87,
    0.0, 0.61, 0.87,
    0.0, 0.61, 0.87,
    1.0,  0.85,  0.0, // left face
    1.0,  0.85,  0.0,
    1.0,  0.85,  0.0,
    1.0,  0.85,  0.0,
    1.0,  0.85,  0.0,
    1.0,  0.85,  0.0,
    0.35,  0.77,  0.0, // top face
    0.35,  0.77,  0.0,
    0.35,  0.77,  0.0,
    0.35,  0.77,  0.0,
    0.35,  0.77,  0.0,
    0.35,  0.77,  0.0,
    0.0,  0.73,  0.49, // bottom face
    0.0,  0.73,  0.49,
    0.0,  0.73,  0.49,
    0.0,  0.73,  0.49,
    0.0,  0.73,  0.49,
    0.0,  0.73,  0.49,
];

var cube_normals = [];

function compute_normals(vertices, normals){
    // Exercise 3: implement the function which computes normals of the cube
    // and add them to cube_normals
}
compute_normals(cube_vertices,cube_normals);

//---------------------------
// definition of the sphere
//---------------------------
var sphere_vertices = [];
var sphere_colors = [];
function create_sphere(){
    let step = 0.01;
    for(let u = 0; u < 1; u = u + step){
        for(let v = 0; v < 1; v = v + step){
            let t = Math.sin(Math.PI*v);

            let x1 = t*Math.cos(2*Math.PI*u);
            let z1 = t*Math.sin(2*Math.PI*u);
            let y1 = Math.cos(Math.PI*v);

            let x4 = t*Math.cos(2*Math.PI*(u+step));
            let z4 = t*Math.sin(2*Math.PI*(u+step));
            let y4 = Math.cos(Math.PI*v);



            t = Math.sin(Math.PI*(v+step));
            let x2 = t*Math.cos(2*Math.PI*u);
            let z2 = t*Math.sin(2*Math.PI*u);
            let y2 = Math.cos(Math.PI*(v+step));

            let x3 = t*Math.cos(2*Math.PI*(u+step));
            let z3 = t*Math.sin(2*Math.PI*(u+step));
            let y3 = Math.cos(Math.PI*(v+step));

            sphere_vertices.push(x1,y1,z1,x3,y3,z3,x2,y2,z2);
            sphere_vertices.push(x1,y1,z1,x4,y4,z4,x3,y3,z3);

            for(let k = 0; k < 6; k++){
                sphere_colors.push(1,1,1);
            }

        }
    }
    //making the sphere a unit sphere
    for(let i = 0; i < sphere_vertices.length; i++){
        sphere_vertices[i] = sphere_vertices[i]/2;
    }
}

create_sphere();

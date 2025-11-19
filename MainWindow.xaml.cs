using System;
using System.IO;
using System.Runtime.InteropServices;
using System.Security.Cryptography.Xml;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using static LineDrawing.MainWindow;

namespace LineDrawing
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    /// 
    public class TreeNode
    {
        public int Value;
        public TreeNode LeftNode;
        public TreeNode RightNode;
    }
    public partial class MainWindow : Window
    {
        [StructLayout(LayoutKind.Sequential, Pack = 1)]
        public struct Pixel
        {
            public byte B;
            public byte G;
            public byte R;

            public static Pixel operator * (Pixel pix, double scalar)
            {
                return new Pixel {B = (byte)(pix.B * scalar), 
                    G = (byte)(pix.G * scalar), 
                    R = (byte)(pix.R * scalar)};
            }

            public static Pixel operator +(Pixel p1, Pixel p2)
            {
                return new Pixel
                {
                    B = (byte)Math.Min(p1.B + p2.B, 255),
                    G = (byte)Math.Min(p1.G + p2.G, 255),
                    R = (byte)Math.Min(p1.R + p2.R, 255)
                };
            }

        }

        public struct Vector3D
        {
            public double X;
            public double Y;
            public double Z;

            public Vector3D(double x, double y, double z)
            {
                X = x;
                Y = y; 
                Z = z;
            }

            public static double DotProduct(Vector3D v1, Vector3D v2)
            {
                return (v1.X * v2.X) + (v1.Y * v2.Y) + (v1.Z * v2.Z);
            }

            public static Vector3D operator -(Vector3D v1, Vector3D v2) {
                return new Vector3D(v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z);
            }

            public static Vector3D operator +(Vector3D a, Vector3D b)
            {
                return new Vector3D(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
            }

            public static Vector3D operator *(Vector3D vector, double scalar)
            {
                return new Vector3D(vector.X * scalar, vector.Y * scalar, vector.Z * scalar);
            }

            public static Vector3D operator *(double scalar, Vector3D vector)
            {
                return vector * scalar; // Reuse the previous overload
            }

            public static Vector3D operator *(Matrix3x3 matrix, Vector3D vector)
            {
                return new Vector3D(matrix.m00 * vector.X + matrix.m01 * vector.Y + matrix.m02 * vector.Z,
                    matrix.m10 * vector.X + matrix.m11 * vector.Y + matrix.m12 * vector.Z,
                    matrix.m20 * vector.X + matrix.m21 * vector.Y + matrix.m22 * vector.Z);
            }

            public static Vector3D operator /(Vector3D vector, double scalar)
            {
                if (scalar == 0)
                {
                    throw new DivideByZeroException("Cannot divide a vector by zero.");
                }
                return new Vector3D(vector.X / scalar, vector.Y / scalar, vector.Z / scalar);
            }

            public double Length()
            {
                return Math.Sqrt(DotProduct(this, this));
            }
        }

        // Matrix3x3 - direct field access
        public struct Matrix3x3
        {
            public double m00, m01, m02;
            public double m10, m11, m12;
            public double m20, m21, m22;

            public Matrix3x3(double v00, double v01, double v02, double v10, double v11, double v12, double v20, double v21, double v22)
            {
                m00 = v00;
                m01 = v01;
                m02 = v02;
                m10 = v10;
                m11 = v11;
                m12 = v12;
                m20 = v20;
                m21 = v21;
                m22 = v22;
            }
        }

        public struct Point
        {
            public int x;
            public int y;

            public Point(int x, int y)
            {
                this.x = x;
                this.y = y;
            }
        }

        public class Scene
        {
            public List<Sphere> spheres;
            public List<Light> lights;
            public Camera camera;
            public double theta;

            public Scene()
            {
                theta = 0;
                spheres = new List<Sphere>();
                lights = new List<Light>();
                camera = new Camera();
            }
        }

        public class Camera
        {
            public Matrix3x3 Rotation;
            public Vector3D Position;
            public Camera()
            {
                Position = new Vector3D(0, 0, 0);
                Rotation = new Matrix3x3
                (
                    1, 0, 0,
                    0, 1, 0,
                    0, 0, 1
                );
            }
        }

        public class Sphere
        {
            public Vector3D Center;
            public double Radius;
            public Pixel Color;
            public double Specular;
            public double Reflective;
        }

        public abstract class Light
        {
            public double Intensity = 0.0;
        }

        public class AmbientLight: Light
        {
        }

        public class PointLight: Light
        {
            public Vector3D Position;
        }

        public class DirectionalLight : Light
        {
            public Vector3D Direction;
        }

        private Scene scene;
        private int Cw = 1000;
        private int Ch = 1000;
        private double Vw = 2.0f;
        private double Vh = 2.0f;
        private double d = 1.0f;
        private static double inf = double.PositiveInfinity;


        public void SaveWriteableBitmapToFile(WriteableBitmap writeableBitmap, string filePath)
        {
            BitmapEncoder encoder;

            // Determine the encoder based on file extension
            string extension = System.IO.Path.GetExtension(filePath).ToLower();
            switch (extension)
            {
                case ".png":
                    encoder = new PngBitmapEncoder();
                    break;
                case ".jpg":
                case ".jpeg":
                    encoder = new JpegBitmapEncoder();
                    break;
                case ".bmp":
                    encoder = new BmpBitmapEncoder();
                    break;
                case ".gif":
                    encoder = new GifBitmapEncoder();
                    break;
                case ".tiff":
                    encoder = new TiffBitmapEncoder();
                    break;
                default:
                    throw new ArgumentException($"Unsupported file extension: {extension}");
            }

            // Create a BitmapFrame from the WriteableBitmap
            BitmapFrame frame = BitmapFrame.Create(writeableBitmap);
            encoder.Frames.Add(frame);

            // Save to file
            using (FileStream stream = new FileStream(filePath, FileMode.Create))
            {
                encoder.Save(stream);
            }
        }

        public MainWindow()
        {
            InitializeComponent();

            scene = new Scene();

            scene.spheres.Add(new Sphere { Center = new Vector3D(0, -1, 3), Radius = 1.0f, 
                Color = new Pixel { R = 255, G = 0, B = 0 }, Specular = 500,
                Reflective = 0.2
            });

            scene.spheres.Add(new Sphere { Center = new Vector3D(2,0,4), Radius = 1.0f, 
                Color = new Pixel { R = 0, G = 0, B = 255 }, Specular = 500,
                Reflective = 0.3
            });

            scene.spheres.Add(new Sphere { Center = new Vector3D(-2,0,4), Radius = 1.0f,
                Color = new Pixel { R = 0, G = 255, B = 0 }, Specular = 10,
                Reflective = 0.4
            });

            scene.spheres.Add(new Sphere { Center = new Vector3D(0, -5001, 0), Radius = 5000.0f,
                Color = new Pixel { R = 255, G = 255, B = 0 }, Specular = 1000,
            Reflective = 0.5});


            scene.lights.Add(new AmbientLight { Intensity = 0.2 });
            scene.lights.Add(new PointLight { Intensity = 0.6, Position = new Vector3D(2, 1, 0) });
            scene.lights.Add(new DirectionalLight { Intensity = 0.2, Direction = new Vector3D(1, 4, 4) });
        }

        private unsafe void DrawLine(int x0, int y0, int x1, int y1, int w, int h, Pixel* buf)
        {
            int dx = Math.Abs(x1 - x0);
            int dy = Math.Abs(y1 - y0);
            int xStep = x0 < x1 ? 1 : -1;
            int yStep = y0 < y1 ? 1 : -1;
            int err = dx - dy;
            var pixel = new Pixel { B = 0, G = 120, R = 30 };
            while (true)
            {
                buf[(y0 * w) + x0] = pixel;
                if (x0 == x1 && y0 == y1)
                    break;
                int err2 = err * 2;
                if (err2 > -dy)
                {
                    err -= dy;
                    x0 += xStep;
                }
                if (err2 < dx)
                {
                    err += dx;
                    y0 += yStep;
                }
            }
        }

        private unsafe void DrawLineEx(int x0, int y0, int x1, int y1, int w, int h, Pixel* buf)
        {
            if (x0 == x1)
            {
                var yStart = Math.Min(y0, y1);
                var yEnd = Math.Max(y0, y1);
                for (var yy = yStart; yy <= yEnd; yy++)
                    buf[yy * w + x0] = new Pixel { B = 0, G = 120, R = 30 };
                return;
            }

            if (y0 == y1)
            {
                var xStart = Math.Min(x0, x1);
                var xEnd = Math.Max(x0, x1);
                for (int x = xStart; x <= xEnd; x++)
                    buf[y0 * w + x] = new Pixel { B = 0, G = 120, R = 30 };
                return;
            }

            // m == slope
            double m = (double)(y1 - y0) / (x1 - x0);
            // y = m.x + b
            double b = y0 - m * x0; 

            if (x0 < x1)
            {
                for (int x = x0; x <= x1; x++)
                {
                    var y = (int)Math.Round(m * x + b);
                    //if (x >= 0 && x < w && y >= 0 && y < h)
                        buf[y * w + x] = new Pixel { B = 0, G = 120, R = 30 };
                }
            }
            else
            {
                for (int x = x1; x >= x0; x--)
                {
                    var y = (int)Math.Round(m * x + b);
                    //if (x >= 0 && x < w && y >= 0 && y < h)
                        buf[y * w + x] = new Pixel { B = 0, G = 120, R = 30 };
                }
            }
        }

        private unsafe double[] Interpolate(int i0, double d0, int i1, double d1) {
            if (i0 == i1) { 
                return new double[] { d0 };
            }
            var values = new List<double>();
            var a = (d1 - d0) / (i1 - i0);
            var d = d0;
            for (int i = i0; i <= i1; i++)
            {
                values.Add(d);
                d = d + a;
            }
            return values.ToArray();
        }

        public unsafe void Swap(Point *p0, Point *p1)
        {
            Point aux = *p0;
            *p0 = *p1;
            *p1 = aux;
        }

        public unsafe void PutPixel(Pixel* image, double x, double y, Pixel color)
        {
            var sx = ((Cw / 2) + (int)Math.Floor(x));
            var sy = ((Ch / 2) - (int)Math.Floor(y)) - 1;

            if (sx == -1)
            {
                sx = 0;
            }
            else if (sx == Cw)
            {
                sx = Cw - 1;
            }

            if (sy == -1)
            {
                sy = 0;
            }
            else if (sy == Ch) 
            { 
                sy = Ch - 1;
            }

            image[(sy * Cw) + sx] = color;
        }

        private unsafe void DrawLineExA(Point p0, Point p1, Pixel color, Pixel* image) {
            if (Math.Abs(p1.x - p0.x) > Math.Abs(p1.y - p0.y)) {
                if (p0.x > p1.x)
                    Swap(&p0, &p1);
                var ys = Interpolate(p0.x, p0.y, p1.x, p1.y);
                for (int x = p0.x; x <= p1.x; x++)
                    PutPixel(image, x, ys[x - p0.x], color);
            }
            else {
                if (p0.y > p1.y)
                    Swap(&p0, &p1);
                var xs = Interpolate(p0.y, p0.x, p1.y, p1.x);
                for (int y = p0.y; y <= p1.y; y++)
                    PutPixel(image, xs[y - p0.y], y, color);
            }
        }

        private unsafe void DrawWireframeTriangle(Point P0, Point P1, Point P2, Pixel color, Pixel* image)
        {
            DrawLineExA(P0, P1, color, image);
            DrawLineExA(P1, P2, color, image);
            DrawLineExA(P2, P0, color, image);
        }

        private unsafe void DrawFilledTriangle(Point P0, Point P1, Point P2, Pixel color, Pixel* image)
        {
            // 1 - Sort the points so that y0 <= y1 <= y2
            if (P1.y < P0.y)
                Swap(&P1, &P0);
            if (P2.y < P0.y)
                Swap(&P2, &P0);
            if (P2.y < P1.y)
                Swap(&P2, &P1);

            double H0 = 0.0f;
            double H1 = 0.0f;
            double H2 = 1.0f;

            // 2 - Compute the x coordinates of the triangle edges
            var X01 = Interpolate(P0.y, P0.x, P1.y, P1.x);
            var H01 = Interpolate(P0.y, H0, P1.y, H1);
            var X12 = Interpolate(P1.y, P1.x, P2.y, P2.x);
            var H12 = Interpolate(P1.y, H1, P2.y, H2);
            var X02 = Interpolate(P0.y, P0.x, P2.y, P2.x);
            var H02 = Interpolate(P0.y, H0, P2.y, H2);

            // 3 - Concatenate the short sides
            X01 = X01.SkipLast(1).ToArray();
            var X012 = X01.Concat(X12).ToArray();
            H01 = H01.SkipLast(1).ToArray();
            var H012 = H01.Concat(H12).ToArray();
            // 4 - Determine which is left and which is right
            var m = (int)Math.Floor((double)(X012.Length / 2));
            double[]? X_LEFT = null;
            double[]? X_RIGHT = null;
            double[]? H_LEFT = null;
            double[]? H_RIGHT = null;
            if (X02[m] < X012[m])
            {
                X_LEFT = X02;
                H_LEFT = H02;
                X_RIGHT = X012;
                H_RIGHT = H012;
            }
            else
            {
                X_LEFT = X012;
                H_LEFT = H012;
                X_RIGHT = X02;
                H_RIGHT = H02;
            }
            // 5 - Draw the horizontal segments
            int y0 = P0.y;
            for (int y = P0.y; y < P2.y; y++)
            {
                var X_L = (int)Math.Ceiling(X_LEFT[y - y0]);
                var X_R = (int)Math.Ceiling(X_RIGHT[y - y0]);
                var H_SEGMENT = Interpolate(X_L, H_LEFT[y-y0], X_R, H_RIGHT[y-y0]);
                for (int x = X_L; x < X_R; x++)
                {
                    var SHADED_COLOR = color * H_SEGMENT[x - X_L];
                    PutPixel(image, x, y, SHADED_COLOR);
                }
            }
        }


        private unsafe Point ViewportToCanvas(double x, double y)
        {
            return new Point((int)Math.Round(x * Cw / Vw) , 
                (int)Math.Round(y * Ch / Vh));
        }

        private unsafe Point ProjectVertex(Vector3D v) {
            return ViewportToCanvas(v.X * d / v.Z, v.Y * d / v.Z);
        }

        private unsafe void DrawSampleCube(Pixel *image)
        {
            // x,y,z
            var vAf = new Vector3D(-1, 1, 2);
            var vBf = new Vector3D(1, 1, 2);
            var vCf = new Vector3D(1, -1, 2);
            var vDf = new Vector3D(-1, -1, 2);

            var vAb = new Vector3D(-1 + 0.35, 1, 2.4);
            var vBb = new Vector3D(1 + 0.35, 1, 2.4);
            var vCb = new Vector3D(1 + 0.35, -1, 2.4);
            var vDb = new Vector3D(-1 + 0.35, -1, 2.4);

            //var vAf = new Vector3D(-1, 1, 2);
            //var vBf = new Vector3D(1, 1, 2);
            //var vCf = new Vector3D(1, -1, 2);
            //var vDf = new Vector3D(-1, -1, 2);

            //var vAb = new Vector3D(-1, 1, 4);
            //var vBb = new Vector3D(1, 1, 4);
            //var vCb = new Vector3D(1, -1, 4);
            //var vDb = new Vector3D(-1, -1, 4);

            // The front face
            var P0 = ProjectVertex(vAf);
            var P1 = ProjectVertex(vBf);
            DrawLineExA(P0, P1, new Pixel { B = 255, G = 10, R = 3 }, image);
            DrawLineExA(ProjectVertex(vBf), ProjectVertex(vCf), new Pixel { B = 255, G = 10, R = 3 }, image);
            DrawLineExA(ProjectVertex(vCf), ProjectVertex(vDf), new Pixel { B = 255, G = 10, R = 3 }, image);
            DrawLineExA(ProjectVertex(vDf), ProjectVertex(vAf), new Pixel { B = 255, G = 10, R = 3 }, image);
            // The back face
            DrawLineExA(ProjectVertex(vAb), ProjectVertex(vBb), new Pixel { B = 3, G = 10, R = 255 }, image);
            DrawLineExA(ProjectVertex(vBb), ProjectVertex(vCb), new Pixel { B = 3, G = 10, R = 255 }, image);
            DrawLineExA(ProjectVertex(vCb), ProjectVertex(vDb), new Pixel { B = 3, G = 10, R = 255 }, image);
            DrawLineExA(ProjectVertex(vDb), ProjectVertex(vAb), new Pixel { B = 3, G = 10, R = 255 }, image);
            // The front-to-back edges
            DrawLineExA(ProjectVertex(vAf), ProjectVertex(vAb), new Pixel { B = 3, G = 255, R = 10 }, image);
            DrawLineExA(ProjectVertex(vBf), ProjectVertex(vBb), new Pixel { B = 3, G = 255, R = 10 }, image);
            DrawLineExA(ProjectVertex(vCf), ProjectVertex(vCb), new Pixel { B = 3, G = 255, R = 10 }, image);
            DrawLineExA(ProjectVertex(vDf), ProjectVertex(vDb), new Pixel { B = 3, G = 255, R = 10 }, image);
        }

        private unsafe void DrawCircleEx(int cx, int cy, int radius, int width, int height, Pixel* img)
        {
            Pixel pixel = new Pixel { B = 0, G = 0, R = 255 };
            for (int i = 1; i <= 360; i++)
            {
                double x = cx - cx;
                double y = (cy - radius) - cy;
                double theta = (1.0f * i) * (Math.PI / 180);
                double xt = (x * Math.Cos(theta)) - (y * Math.Sin(theta));
                double yt = (x * Math.Sin(theta)) + (y * Math.Cos(theta));
                xt += cx;
                yt += cy;
                x = Math.Round(xt);
                y = Math.Round(yt);
                img[(int)((y * width) + x)] = pixel;
            }
        }

        private unsafe void RenderClick(object sender, RoutedEventArgs e)
        {
            var wbmp = new WriteableBitmap(Cw, Ch, 96.0f, 96.0f, PixelFormats.Bgr24, null);
            var mem = (Pixel*)wbmp.BackBuffer;
            wbmp.Lock();
            var pixel = new Pixel { B = 255, G = 255, R = 255 };
            for (int y = 0; y < Ch; y++)
                for (int x = 0; x < Cw; x++)
                    mem[y * Cw + x] = pixel;
            //DrawLineExA(new Point(50, 50), new Point(81, 75), w, h, new Pixel { B = 33, G = 25, R = 4 },mem);
            //DrawLineExA(new Point(250, 300), new Point(800, 750), w, h, new Pixel { B = 33, G = 25, R = 4 }, mem);
            //DrawLineExA(new Point(81, 75), new Point(50, 700), w, h, new Pixel { B = 33, G = 25, R = 4 }, mem);

            var P0 = new Point(-200, -250);
            var P1 = new Point(200, 50);
            var P2 = new Point(20, 250);
            var color = new Pixel { R = 255, G = 10, B = 3 };
            var fillColor = new Pixel { R = 3, G = 190, B = 3 };
            //DrawWireframeTriangle(P0, P1, P2, color, mem);
            //DrawFilledTriangle(P0, P1, P2, fillColor, mem);
            DrawSampleCube(mem);
            wbmp.AddDirtyRect(new Int32Rect(0, 0, Cw, Ch));
            wbmp.Unlock();
            img.Source = wbmp;
        }

        private unsafe void GraphicsClick(object sender, RoutedEventArgs e)
        {
            var w = 1920;
            var h = 1080;
            
            var wbmp = new WriteableBitmap(w, h, 96.0f, 96.0f, PixelFormats.Bgr24, null);
            var mem = (Pixel*)wbmp.BackBuffer;            
            wbmp.Lock();

            var i = 1;
            var hp = h / 4;
            var wp = w / 4;
            var rng = new Random();
            for (int hc = 0; hc < 4; hc++)
                for (int wc = 0; wc < 4; wc++)
                {
                    var pixel = new Pixel { R = (byte)(i * 4), G = 0, B = (byte)(i * 8) };
                    for (int y = (hp * (hc + 0)); y < (hp * (hc + 1)); y++)
                        for (int x = (wp * (wc + 0)); x < (wp * (wc + 1)); x++)
                            mem[w * y + x] = pixel;
                    i++;
                }

            wbmp.AddDirtyRect(new Int32Rect(0, 0, w, h));
            wbmp.Unlock();

            img.Source = wbmp;
        }

        private unsafe Vector3D CanvasToViewport(double x, double y) {
            return new Vector3D(x*Vw/Cw, y*Vh/Ch, d);
        }

        private Tuple<double, double> IntersectRaySphere(Vector3D O, Vector3D D, Sphere sphere)
        {
            var r = sphere.Radius;
            Vector3D CO = O - sphere.Center;

            var a = Vector3D.DotProduct(D, D);
            var b = 2 * Vector3D.DotProduct(CO, D);
            var c = Vector3D.DotProduct(CO, CO) - r * r;

            var discriminant = b * b - 4 * a * c;
            if (discriminant < 0) {
                return new Tuple<double, double>(inf, inf);
            }

            var t1 = (-b + Math.Sqrt(discriminant)) / (2 * a);
            var t2 = (-b - Math.Sqrt(discriminant)) / (2 * a);

            return new Tuple<double, double>(t1, t2);
        }

        private double ComputeLighting(Vector3D P, Vector3D N, Vector3D V, double s)
        {
            double i = 0.0;
            Vector3D L;
            double t_max;

            foreach (var light in scene.lights) 
            {
                if (light is null)
                    continue;
                if (light is AmbientLight)
                {
                    i += ((AmbientLight)light).Intensity;
                }
                else
                {                    
                    if (light is PointLight)
                    {
                        L = (light as PointLight).Position - P;
                        t_max = 1.0;
                    }
                    else
                    {
                        L = (light as DirectionalLight).Direction;
                        t_max = inf;
                    }

                    var tup = ClosestIntersection(P, L, 0.001, t_max);
                    Sphere shadow_sphere = tup.Item1;
                    double shadow_t = tup.Item2;

                    if (shadow_sphere != null) {
                        continue;
                    }

                    var n_dot_l = Vector3D.DotProduct(N, L);
                    if (n_dot_l > 0)
                    {
                        i += light.Intensity * n_dot_l / (N.Length() * L.Length());
                    }

                    if (s != -1.00) {
                        var R = (2 * N) * n_dot_l - L;
                        var r_dot_v = Vector3D.DotProduct(R, V);
                        if (r_dot_v > 0) {
                            i += light.Intensity * Math.Pow(r_dot_v / (R.Length() * V.Length()), s);
                        }
                    }
                }
            }
            if (i > 1)
                i = 1;
            return i;
        }

        private unsafe Tuple<Sphere, double> ClosestIntersection(Vector3D O, Vector3D D, double t_min, double t_max)
        {
            double closest_t = inf;
            Sphere closest_sphere = null;


            foreach (var sphere in scene.spheres)
            {
                Tuple<double, double> ret = IntersectRaySphere(O, D, sphere);

                double t1 = ret.Item1;
                double t2 = ret.Item2;

                if ((t1 >= t_min && t1 <= t_max) && t1 < closest_t)
                {
                    closest_t = t1;
                    closest_sphere = sphere;
                }

                if ((t2 >= t_min && t2 <= t_max) && t2 < closest_t)
                {
                    closest_t = t2;
                    closest_sphere = sphere;
                }
            }
            return new Tuple<Sphere, double>(closest_sphere, closest_t);
        }

        private unsafe Vector3D ReflectRay(Vector3D R, Vector3D N) {
            return 2 * N * Vector3D.DotProduct(N , R) - R;
        }

        private unsafe Pixel TraceRay(Vector3D O, Vector3D D, double t_min, double t_max, int recursion_depth) {
            var closest = ClosestIntersection(O, D, t_min, t_max);

            Sphere closest_sphere = closest.Item1;
            double closest_t = closest.Item2;

            if (closest_sphere == null)
            {
                return new Pixel { R = 255, G = 255, B = 255 };
            }

            var P = O + closest_t * D;
            var N = P - closest_sphere.Center;
            N = N / N.Length();
            var local_color = closest_sphere.Color * ComputeLighting(P, N, -1 * D, closest_sphere.Specular);
            
            var r = closest_sphere.Reflective;
            if (recursion_depth <= 0 || r <= 0) {
                return local_color;
            }
            var R = ReflectRay(-1 * D, N);
            var reflected_color = TraceRay(P, R, 0.001, inf, recursion_depth - 1);
            return (local_color * (1 - r)) + reflected_color * r;
        }

        private unsafe void RayTraceClick(object sender, RoutedEventArgs e)
        {

            var wbmp = new WriteableBitmap(Cw, Ch, 96.0f, 96.0f, PixelFormats.Bgr24, null);
            var mem = (Pixel*)wbmp.BackBuffer;
            wbmp.Lock();
            double cos = Math.Cos(Math.PI / 180 * scene.theta);
            double sin = Math.Sin(Math.PI / 180 * scene.theta);
            scene.camera.Rotation = new Matrix3x3
            (
                cos, 0, sin,
                0, 1, 0,
                -sin, 0, cos
            );

            Vector3D O = new Vector3D(0, 0, 0);
            for (int x = -Cw / 2; x < Cw / 2; x++)
            {
                for (int y = -Ch / 2; y < Ch / 2; y++)
                {
                    var D = scene.camera.Rotation * CanvasToViewport(x, y);
                    Pixel color = TraceRay(O, D, 1, inf, 3);
                    PutPixel(x, y, color, mem);
                }
            }
            wbmp.AddDirtyRect(new Int32Rect(0, 0, Cw, Ch));
            wbmp.Unlock();
            img.Source = wbmp;
            scene.theta += 1.0;
            //SaveWriteableBitmapToFile(wbmp, "E:\\Dump\\raytrace.bmp");
        }

        private unsafe void PutPixel(int x, int y, Pixel pixel, Pixel* mem)
        {
            var sx = ((Cw / 2) + x);
            var sy = ((Ch / 2) - y) - 1;
            mem[(sy * Cw) + sx] = pixel;
        }
    }
}
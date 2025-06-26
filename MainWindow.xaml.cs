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

        public class Scene
        {
            public List<Sphere> spheres;
            public List<Light> lights;

            public Scene()
            {
                spheres = new List<Sphere>();
                lights = new List<Light>();
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
        private double Vw = 1.0f;
        private double Vh = 1.0f;
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
            var w = 320;
            var h = 240;
            var wbmp = new WriteableBitmap(w, h, 96.0f, 96.0f, PixelFormats.Bgr24, null);
            var mem = (Pixel*)wbmp.BackBuffer;
            wbmp.Lock();
            var pixel = new Pixel { B = 255, G = 100, R = 48 };
            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                    mem[y * w + x] = pixel;
            //DrawLineEx(50, 50, 100, 100, w, h, mem);
            //DrawLineEx(50, 50, 200, 200, w, h, mem);
            //DrawLineEx(75, 75, 100, 50, w, h, mem);
            //DrawLineEx(100, 50, 125, 75, w, h, mem);
            //DrawLineEx(125, 75, 75, 75, w, h, mem);
            DrawLineEx(50, 50, 81, 75, w, h, mem);
            //DrawCircleEx(100, 100, 30, w, h, mem);
            wbmp.AddDirtyRect(new Int32Rect(0, 0, w, h));
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

        private unsafe Pixel TraceRay(Vector3D O, Vector3D D, double t_min, double t_max) {
            var closest = ClosestIntersection(O, D, t_min, t_max);

            Sphere closest_sphere = closest.Item1;
            double closest_t = closest.Item2;

            if (closest_sphere == null)
            {
                return new Pixel { R = 255, G = 255, B = 255 };
            }

            var P = O + closest_t * D;
            var N = P - closest_sphere.Center;
            var len = N.Length();
            N = N / len;
            var lighting = ComputeLighting(P, N, -1 * D, closest_sphere.Specular);
            var colorAppliedLighting = closest_sphere.Color * lighting;
            return colorAppliedLighting;
        }

        private unsafe void RayTraceClick(object sender, RoutedEventArgs e)
        {
            var wbmp = new WriteableBitmap(Cw, Ch, 96.0f, 96.0f, PixelFormats.Bgr24, null);
            var mem = (Pixel*)wbmp.BackBuffer;
            wbmp.Lock();

            Vector3D O = new Vector3D(0, 0, 0);
            for (int x = -Cw / 2; x < Cw / 2; x++) {
                for (int y = -Ch / 2; y < Ch / 2; y++)
                {
                    var D = CanvasToViewport(x, y);
                    Pixel color = TraceRay(O, D, 1, inf);
                    PutPixel(x, y, color, mem);
                }
            }

            wbmp.AddDirtyRect(new Int32Rect(0, 0, Cw, Ch));
            wbmp.Unlock();

            img.Source = wbmp;
            SaveWriteableBitmapToFile(wbmp, "E:\\Dump\\raytrace.bmp");
        }

        private unsafe void PutPixel(int x, int y, Pixel pixel, Pixel* mem)
        {
            // X >= -960 && X <= 960
            // Y >= -540 && Y <= 540
            var sx = ((Cw / 2) + x);
            var sy = ((Ch / 2) - y) - 1;
            mem[(sy * Cw) + sx] = pixel;
            //0 -> 6 220 799
        }
    }
}
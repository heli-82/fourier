use num::complex::Complex64;
use std::{f64::consts::PI, vec};
use plotly::box_plot::{BoxMean, BoxPoints};
use plotly::common::{ErrorData, ErrorType, Line, Marker, Mode, Orientation, Title};
use plotly::histogram::{self, Bins, Cumulative, HistFunc, HistNorm};
use plotly::layout::{Axis, BarMode, BoxMode, Layout, Margin};
use plotly::{Bar, BoxPlot, Histogram, Plot, color::{NamedColor, Rgb, Rgba}, Scatter};
fn dft(input: &[f64]) -> Vec<Complex64> {
    let n = input.len();
    println!("{:?}", n);

    let mut out: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n];

    for k in 0..n {
        let mut sum = Complex64::new(0.0, 0.0);
        for l in 0..n {
            let angle: f64 = (-2.0) * PI * (k as f64) * (l as f64) / (n as f64);
            let wn = Complex64::new(
                if angle.cos().abs() < 1e-10 {
                    0.0
                } else {
                    angle.cos()
                },
                if angle.sin().abs() < 1e-10 {
                    0.0
                } else {
                    angle.sin()
                },
            );
            sum += input[l] * wn;
            sum = Complex64::new(
                if sum.re.abs() > 1e-10 { sum.re } else { 0.0 },
                if sum.im.abs() > 1e-10 { sum.im } else { 0.0 },
            );
        }
        out[k] = sum;
    }
    out
}

fn idft(spectrum: &[Complex64]) -> Vec<Complex64> {
    let n = spectrum.len();
    let mut out: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n];

    for l in 0..n {
        let mut sum = Complex64::new(0.0, 0.0);
        for k in 0..n {
            let angle: f64 = 2.0 * PI * (k as f64) * (l as f64) / (n as f64);
            let wn = Complex64::new(angle.cos(), angle.sin());
            sum += spectrum[k] * wn;
        }
        out[l] = sum / (n as f64);
    }
    out
}

fn fft(signal: &mut [Complex64]) {
    let n = signal.len();

    if n <= 0 {
        return;
    }

    let (mut even, mut odd) = (
        vec![Complex64::new(0.0, 0.0); n / 2],
        vec![Complex64::new(0.0, 0.0); n / 2],
    );
    for i in 0..(n / 2) {
        even[i] = signal[2 * i];
        odd[i] = signal[2 * i + 1];
    }

    fft(&mut even);
    fft(&mut odd);

    for k in 0..((n as usize) / 2) {
        let mut t = Complex64::from_polar(1.0, (-2.0) * PI * (k as f64) / (n as f64)) * odd[k];
        t = Complex64::new(
            if t.re.abs() > 1e-10 { t.re } else { 0.0 },
            if t.im.abs() > 1e-10 { t.im } else { 0.0 },
        );
        signal[k] = even[k] + t;
        signal[k + n / 2] = even[k] - t;
    }
}

fn main() {
    let a = vec![0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0];
    println!("{:?}", dft(&a));

    /*let mut b: Vec<Complex64> = vec![
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
    ];
    fft(&mut b);
    println!("{:?}", b);*/

    let mut c: Vec<Complex64> = a.iter().map(|x| Complex64::new(*x, 0.0)).collect();
    fft(&mut c);
    println!("{:?}", c);


    // let samples = 
    
    // let trace = Histogram::new(
    //     (0..(c.len())).collect::<Vec<usize>>(),
    //     // c.iter().map(|x| (x.re*x.re + x.im*x.im).sqrt()).collect::<Vec<f64>>(),
    // );
    let trace = histogram::Histogram::new_xy((0..(c.len())).map(|x| x as f64).collect::<Vec<f64>>(), c.iter().map(|x| (x.re*x.re + x.im*x.im).sqrt()).collect::<Vec<f64>>());
    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.show();
}

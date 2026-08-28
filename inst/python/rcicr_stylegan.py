"""Render latents through a StyleGAN for rcicr's latentGeneratorCommand().

rcicr calls this script as

    python rcicr_stylegan.py --network <pkl> [options] <latents.csv> <outdir>

reads back the PNGs it writes, and never loads the model itself. That is the
whole interface: anything able to honour it can stand in for this script.

The latents file has one comma-separated row per stimulus and no header. One
PNG per row goes into the output directory, named 00001.png upwards in the
order of the rows.

This needs a StyleGAN2 or StyleGAN3 checkout on sys.path, plus torch. Neither
is an rcicr dependency, because rcicr only ever runs the command.

Point rcicr at it with, in R:

    generator <- latentGeneratorCommand(
      command    = "python",
      args       = c(system.file("python", "rcicr_stylegan.py", package = "rcicr"),
                     "--network", "ffhq.pkl", "--space", "w"),
      latent_dim = 512,
      img_size   = 1024,
      space      = "w",
      latent_mean = <the average W latent>,
      latent_sd   = <the spread of W along each dimension>,
      weights     = "ffhq.pkl"
    )

latent_mean and latent_sd are what perturbation sizes are measured against, so
they decide what latent_sigma = 1 means. --report-stats prints values for them,
sampled from the network, in a form R can paste in.
"""

import argparse
import csv
import pickle
import sys

import numpy as np


def load_network(path, device):
    import torch  # noqa: F401  (imported for its side effect on the pickle)

    with open(path, "rb") as handle:
        data = pickle.load(handle)

    return data["G_ema"].to(device).eval()


def read_latents(path):
    with open(path, newline="") as handle:
        rows = [[float(value) for value in row] for row in csv.reader(handle) if row]

    return np.asarray(rows, dtype=np.float32)


def to_uint8(image):
    """StyleGAN emits roughly [-1, 1]; PNG wants [0, 255]."""
    return (image.clamp(-1, 1) * 127.5 + 128).clamp(0, 255).to("cpu", dtype=int)


def render(network, latents, space, device, truncation):
    import torch

    batch = torch.from_numpy(latents).to(device)

    with torch.no_grad():
        if space == "z":
            label = torch.zeros([batch.shape[0], network.c_dim], device=device)
            return network(batch, label, truncation_psi=truncation, noise_mode="const")

        # W is one vector per image; the synthesis network wants it repeated
        # once per layer, which is what W+ already is.
        if space == "w":
            batch = batch.unsqueeze(1).repeat([1, network.num_ws, 1])

        return network.synthesis(batch, noise_mode="const")


def report_stats(network, device, samples, space):
    import torch

    with torch.no_grad():
        z = torch.randn([samples, network.z_dim], device=device)
        label = torch.zeros([samples, network.c_dim], device=device)
        w = network.mapping(z, label)[:, 0, :] if space != "z" else z

    w = w.cpu().numpy()
    fmt = lambda values: ", ".join("%.6f" % v for v in values)  # noqa: E731

    print("latent_mean <- c(%s)" % fmt(w.mean(axis=0)))
    print("latent_sd <- c(%s)" % fmt(w.std(axis=0)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", required=True, help="path to the .pkl checkpoint")
    parser.add_argument("--space", default="w", choices=["z", "w", "w+"])
    parser.add_argument("--truncation", type=float, default=1.0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--batch", type=int, default=8)
    parser.add_argument("--report-stats", type=int, default=0, metavar="N",
                        help="print latent_mean and latent_sd from N samples, then exit")
    parser.add_argument("latents", nargs="?", help="the latents CSV rcicr wrote")
    parser.add_argument("outdir", nargs="?", help="where to write the PNGs")
    args = parser.parse_args(argv)

    import torch
    from PIL import Image

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    network = load_network(args.network, device)

    if args.report_stats:
        report_stats(network, device, args.report_stats, args.space)
        return 0

    if not args.latents or not args.outdir:
        parser.error("latents and outdir are required unless --report-stats is given")

    latents = read_latents(args.latents)

    written = 0
    for start in range(0, len(latents), args.batch):
        images = render(network, latents[start:start + args.batch], args.space,
                        device, args.truncation)
        for image in images:
            written += 1
            array = to_uint8(image.permute(1, 2, 0)).numpy().astype("uint8")
            Image.fromarray(array).save("%s/%05d.png" % (args.outdir, written))

    return 0


if __name__ == "__main__":
    sys.exit(main())

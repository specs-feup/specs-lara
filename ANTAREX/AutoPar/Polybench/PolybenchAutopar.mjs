import Query from "@specs-feup/lara/api/weaver/Query.js";
import { FileJp, FunctionJp, Loop } from "@specs-feup/clava/api/Joinpoints.js";
import Parallelize from "@specs-feup/clava-autopar/api/Parallelize.js";
import AutoParStats from "@specs-feup/clava-autopar/api/AutoParStats.js";

// Reset stats
AutoParStats.reset();

const $loops = Query.search(FileJp)
    .search(FunctionJp, {
        // Only parallelize loops inside Polybench kernel functions
        name: (name) => name.startsWith("kernel_"),
    })
    .search(Loop)
    .get();

// Set name
$loops.forEach(($loop) =>
    AutoParStats.get().setName(
        Io.removeExtension($loop.getAncestor("file").name)
    )
);

Parallelize.forLoops($loops);

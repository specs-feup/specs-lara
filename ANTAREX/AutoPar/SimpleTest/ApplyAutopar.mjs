import Query from "@specs-feup/lara/api/weaver/Query.js"
import { FileJp, Loop } from "@specs-feup/clava/api/Joinpoints.js"
import Parallelize from "@specs-feup/clava-autopar/api/Parallelize.js"

const $loops = Query.search(FileJp).search(Loop).get();

Parallelize.forLoops($loops);

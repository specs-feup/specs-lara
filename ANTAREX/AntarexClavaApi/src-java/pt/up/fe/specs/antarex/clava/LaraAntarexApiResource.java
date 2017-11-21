/**
 * Copyright 2013 SuikaSoft.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
 * an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License. under the License.
 */

package pt.up.fe.specs.antarex.clava;

import org.lara.interpreter.weaver.utils.LaraResourceProvider;

/**
 * @author Joao Bispo
 *
 */
public enum LaraAntarexApiResource implements LaraResourceProvider {
    TEST("Test.lara"),

    // MultiVersioning
    MULTI_POINTERS("multi/MultiVersionPointers.lara"),
    MULTI_POINTERS_ASPECTS("multi/MultiVersionPointersAspects.lara"),

    // LIBVC
    LIBVC("libvc/LibVC.lara"),
    LIBVC_ASPECTS("libvc/LibVCAspects.lara"),

    // MARGOT
    MARGOT_CODE_GEN("margot/MargotCodeGen.lara"),
    MARGOT_CODE_GEN_ASPECTS("margot/MargotCodeGenAspects.lara"),
    MARGOT_CODE_GEN_STRINGS("margot/MargotStringsGen.lara"),
    // MargotConfig and all the sub files that are imported by MargotConfig
    MARGOT_CONFIG("margot/MargotConfig.lara"),
    MARGOT_CONFIG_BLOCK("margot/MargotConfigImports/MargotBlock.lara"),
    MARGOT_CONFIG_MONITOR("margot/MargotConfigImports/MargotMonitor.lara"),
    MARGOT_CONFIG_ENERGY_MONITOR("margot/MargotConfigImports/MargotEnergyMonitor.lara"),
    MARGOT_CONFIG_TIME_MONITOR("margot/MargotConfigImports/MargotTimeMonitor.lara"),
    MARGOT_CONFIG_THROUGHPUT_MONITOR("margot/MargotConfigImports/MargotThroughputMonitor.lara"),
    MARGOT_CONFIG_STATE("margot/MargotConfigImports/MargotState.lara");

    private final String resource;

    private static final String WEAVER_PACKAGE = "clava/";
    // This is the prefix that will appear in the LARA import, WEAVER_PACKAGE will be ignored
    // E.g., import antarex.Test;
    private static final String BASE_PACKAGE = "antarex/";

    /**
     * @param resource
     */
    private LaraAntarexApiResource(String resource) {
        this.resource = WEAVER_PACKAGE + getSeparatorChar() + BASE_PACKAGE + resource;
    }

    /* (non-Javadoc)
     * @see org.suikasoft.SharedLibrary.Interfaces.ResourceProvider#getResource()
     */
    @Override
    public String getOriginalResource() {
        return resource;
    }

}

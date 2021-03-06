#pragma once

#define STR_EXPAND(tok) # tok
#define STR(tok) STR_EXPAND(tok)
#define PAW_VERSION_MAJOR @PAW_VERSION_MAJOR@
#define PAW_VERSION_MINOR @PAW_VERSION_MINOR@
#define PAW_VERSION_PATCH @PAW_VERSION_PATCH@
#define PAW_GIT_BRANCH @PAW_GIT_BRANCH@
#define PAW_GIT_COMMIT_SHORT_HASH @PAW_GIT_COMMIT_SHORT_HASH@
#define PAW_GIT_COMMIT_LONG_HASH @PAW_GIT_COMMIT_LONG_HASH@
#define PAW_GIT_NUM_DIRTY_FILES @PAW_GIT_NUM_DIRTY_FILES@
#define PAW_SOURCE_DIR @PROJECT_SOURCE_DIR@
#define PAW_BINARY_DIR @PROJECT_BINARY_DIR@
#define PAW_BOOST_FOUND @Boost_FOUND@

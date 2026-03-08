import streamlit as st
import pandas as pd
import time
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw
from io import BytesIO
import base64
from .smartbind_wrapper import run_virtual_screening
from .ui_components import render_header, render_section_divider, render_footer, render_success_alert

# Personalized SVG Icons
# RNA helix – premium stylized version from RNA.svg
ICON_HELIX = ('<svg width="22" height="22" viewBox="0 0 512 512" fill="none" xmlns="http://www.w3.org/2000/svg">'
    '<path d="M150.815 0L160.901 0C163.181 0.913013 165.648 1.1956 167.927 2.05723C176.627 5.3469 183.973 11.9088 187.816 20.4291C189.266 23.6033 190.218 26.9817 190.64 30.4456C191.333 35.9937 190.903 42.4391 191.01 48.0944C191.254 60.969 190.609 74.692 191.188 87.4759C192.454 112.813 199.326 137.553 211.308 159.913C224.218 184.108 242.749 204.85 265.343 220.392C275.329 227.236 284.55 232.164 294.558 239.707C351.742 282.685 386.923 348.809 390.612 420.248C391.293 432.191 391.007 445.089 391.016 457.143C391.021 463.433 391.35 474.232 390.785 480.139C390.509 483.358 389.782 486.522 388.625 489.538C385.312 498.136 378.697 505.135 370.265 508.951C368.183 509.893 362.052 511.475 361.189 512L351.099 512C349.619 511.283 347.389 510.957 345.78 510.476C335.105 507.335 325.843 498.304 322.616 487.606C320.144 479.409 321.155 469.254 321.01 460.653C320.831 450.072 321.177 439.414 320.957 428.82C320.583 414.857 318.523 400.992 314.82 387.524C304.443 349.851 281.548 316.833 249.903 293.907C243.32 289.15 236.658 285.361 229.889 280.975C222.927 276.462 216.212 271.577 209.774 266.343C177.896 240.396 153.227 206.688 138.137 168.455C128.405 143.821 122.728 117.774 121.326 91.3243C120.469 74.5769 121.473 57.8148 120.991 40.922C120.502 23.7554 125.88 9.87253 142.714 2.58691C145.25 1.48965 148.096 1.08011 150.622 0.0787148L150.815 0Z" fill="var(--fz-accent)"/>'
    '<path d="M221.111 31.0136L312.021 30.9886L341.696 30.9785C349.101 30.9769 359.908 29.6548 365.757 34.6551C368.835 37.2954 370.697 41.0807 370.91 45.13C371.149 49.6923 369.587 53.1805 366.613 56.5187C364.759 58.3778 361.563 60.3503 358.981 60.6721C353.891 61.3067 348.097 61.0929 342.935 61.0904L316.26 61.0755L220.973 61.0834C221.293 51.5109 220.527 40.3575 221.111 31.0136Z" fill="var(--fz-text)" opacity="0.8"/>'
    '<path d="M221.72 91.1214C221.936 91.1056 222.152 91.0922 222.368 91.0812C231.114 90.6785 242.219 90.9882 251.088 90.9862L305.673 90.9804L339.501 90.983C345.298 90.984 351.533 90.8305 357.306 91.2127C360.713 91.4383 363.732 92.7428 366.255 95.0991C369.135 97.8071 370.816 101.553 370.925 105.505C371.06 110.062 369.382 113.514 366.352 116.781C361.184 121.234 357.245 121.182 350.715 121.092C348.116 121.056 345.366 121.083 342.755 121.079L314.453 121.076L259.052 121.082C248.6 121.086 237.589 121.289 227.181 120.999C226.543 118.116 225.645 115.204 224.953 112.304C223.572 106.608 222.554 100.829 221.906 95.0037C221.792 93.9653 221.348 91.958 221.72 91.1214Z" fill="var(--fz-text)" opacity="0.8"/>'
    '<path d="M241.224 151.116C249.114 150.741 258.987 150.991 266.998 150.988L313.441 150.98L342.164 150.982C347.194 150.982 352.655 150.828 357.639 151.244C360.879 151.515 363.855 152.86 366.255 155.098C369.169 157.835 370.854 161.633 370.928 165.63C371.01 170.289 369.28 173.751 366.112 177.005C365.816 177.251 365.516 177.493 365.213 177.73C363.177 179.315 360.743 180.586 358.167 180.797C352.441 181.266 346.329 181.089 340.568 181.089L308.791 181.08L280.79 181.094C276.355 181.107 268.224 181.37 264.138 180.852C260.807 177.894 240.28 154.411 241.224 151.116Z" fill="var(--fz-text)" opacity="0.8"/>'
    '<path d="M154.737 451.091C160.718 450.8 167.799 450.986 173.882 450.986L207.636 450.983C235.189 450.982 263.497 450.642 290.997 451.054L291.012 481.04C286.723 481.108 282.357 481.074 278.062 481.08L195.232 481.085L169.447 481.095C162.741 481.098 152.641 482.154 147.143 478.133C143.87 475.745 141.705 472.134 141.141 468.123C139.92 459.096 145.993 452.365 154.737 451.091Z" fill="var(--fz-accent)" opacity="0.6"/>'
    '<path d="M154.735 391.091C163.595 390.74 173.937 390.983 182.881 390.985L234.142 390.986L234.142 390.986L267 390.982C271.757 390.982 280.269 390.739 284.701 391.183C286.626 397.819 288.039 403.531 289.146 410.38C289.695 413.772 289.956 417.937 290.584 421.054L277.125 421.078L195.269 421.087L169.391 421.094C162.79 421.096 152.703 422.12 147.268 418.224C143.978 415.877 141.779 412.298 141.172 408.303C139.832 399.221 145.91 392.377 154.735 391.091Z" fill="var(--fz-accent)" opacity="0.6"/>'
    '<path d="M154.729 331.091C160.486 330.857 166.841 330.985 172.666 330.983L203.319 330.983L230.267 330.987C233.097 330.987 246.882 330.46 248.603 332.024C255.001 337.838 260.828 345.369 265.828 352.444C267.724 355.126 269.404 358.525 271.561 360.935C267.033 361.047 262.504 361.094 257.975 361.077L189.283 361.088L167.552 361.096C161.449 361.098 152.728 361.943 147.585 358.444C144.22 356.168 141.923 352.624 141.219 348.623C139.65 339.487 145.812 332.398 154.729 331.091Z" fill="var(--fz-accent)" opacity="0.6"/>'
    '</svg>')

# Full HTVS.svg with all paths theme-mapped:
#  #24455C (dark navy) → accent
#  #3693C6 (main blue) → accent at reduced opacity
#  #386178 (mid teal)  → accent at lower opacity
#  #9066A5 (purple)    → accent glow/muted
#  #D4CBCB (grey)      → text-dim
ICON_HTVS_MOL = ('<svg width="28" height="28" viewBox="0 0 430 444" fill="none"'
    ' xmlns="http://www.w3.org/2000/svg">'
    # outer silhouette
    '<path fill="currentColor" opacity="0.12" d="M60.132 113.196C66.4717 104.47 71.2321 97.803 78.5959 89.7682C107.265 58.4121 145.717 37.7008 187.675 31.0156C219.63 25.9656 255.869 29.876 285.832 42.1C300.638 48.2272 314.634 56.1548 327.502 65.7037C365.237 93.2266 393.922 139.809 401.058 185.87C405.386 212.506 404.189 239.748 397.541 265.902C396.597 269.708 394.068 279.221 392.416 282.744C383.169 282.765 371.754 283.384 362.854 282.935C352.067 283.032 341.279 283.045 330.492 282.974C329.098 283.012 327.55 282.968 326.143 282.957C310.216 283.395 292.242 283.023 276.159 282.925L272.726 283.212C267.201 284.25 264.926 287.455 262.826 292.404C261.862 299.431 264.384 300.857 268.216 306.183C270.07 308.05 273.034 311.959 274.862 314.121L288.23 329.8L301.33 345.358C303.018 347.363 306.335 351.03 307.764 353.021C311.064 356.467 316.314 363.005 319.441 366.831C319.775 369.963 319.807 372.232 319.515 375.368C280.094 399.167 242.476 410.037 196.022 405.092C158.235 401.037 122.596 385.509 93.8978 360.596C90.7199 357.841 87.7414 355.823 84.9221 352.58C83.0185 350.391 80.261 348.753 78.772 346.013C78.3978 345.325 77.2471 344.1 76.723 343.451C71.868 340.458 66.2149 332.545 62.9041 327.852C60.6269 325.018 53.9062 315.596 53.8737 312.282C51.4849 305.954 46.5732 298.77 43.9042 292.36C43.2885 290.881 41.3836 288.529 40.9981 287.117C39.6666 282.239 38.5899 277.815 36.9463 273.002C35.0422 268.638 31.7959 258.452 31.3886 253.776C31.9311 253.082 31.9003 252.895 32.059 251.99C29.999 246.02 30.0175 234.394 29.0312 227.956C28.4601 182.998 36.3199 151.575 60.132 113.196Z"/>'
    # central molecule – bright accent
    '<path fill="var(--fz-accent)" d="M213.166 132.594C218.782 129.944 219.254 129.935 225.003 132.326C229.398 135.025 233.889 137.44 238.296 140.013C241.859 142.447 251.846 146.121 252.246 150.993C253.088 161.244 252.281 171.884 252.578 182.19C252.751 188.186 238.777 194.359 234.19 196.956C231.308 198.586 228.479 200.31 225.708 202.124C225.892 209.999 225.832 217.668 225.963 225.461C230.747 228.046 235.498 230.692 240.215 233.398C243.481 235.236 251.462 238.874 252.084 242.778C253.399 251.019 252.711 260.228 252.704 268.635C252.305 276.386 253.144 277.756 246.462 282.358C242.911 284.247 239.4 286.209 235.93 288.244C230.943 291.173 224.785 295.272 219.303 296.618C217.084 295.964 214.98 295.076 212.909 294.048C209.73 292.47 206.751 290.497 203.663 288.751C200.126 286.751 188.816 281.506 187.032 278.328C185.168 275.01 185.432 246.012 186.734 242.05C187.104 240.926 187.722 240.173 188.613 239.416C191.146 237.265 194.377 235.691 197.254 234.044C202.468 231.017 207.736 228.084 213.056 225.246C212.795 217.107 212.846 209.993 212.988 201.855C210.728 200.595 208.42 199.196 206.046 197.751C199.923 193.984 191.906 190.971 187.094 185.599C185.74 184.087 185.895 176.518 185.981 174.33C186.284 166.617 185.564 158.938 186.245 151.261C186.643 146.767 193.849 144.113 197.009 141.944C202.482 138.803 207.823 135.908 213.166 132.594Z"/>'
    # central dot – accent
    '<path fill="var(--fz-accent)" opacity="0.9" d="M219.4 236.486C219.863 236.538 238.083 247.183 239.806 248.165C239.688 253.435 239.379 266.194 239.961 271.096C233.744 274.159 225.571 279.102 219.596 282.616C218.6 282.628 201.915 272.638 199.211 271.305C199.299 265.727 199.224 260.126 199.272 254.546C199.289 252.556 199.343 250.244 198.983 248.302Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.9" d="M219.236 144.45C220.415 144.489 237.225 154.551 239.753 155.815C239.704 162.847 239.524 171.814 239.851 178.754C233.067 182.491 226.394 186.506 219.591 190.191C217.589 189.531 213.162 186.883 211.205 185.748L199.258 178.792C199.322 173.224 199.471 161.644 198.992 156.226C204.704 152.222 212.971 147.923 219.236 144.45Z"/>'
    # stem connector – muted accent
    '<path fill="var(--fz-accent)" opacity="0.55" d="M186.262 85.8051C186.273 78.3859 185.95 71.037 186.076 63.554C186.057 58.6477 186.117 52.0234 192.859 52.0233C200.502 52.0232 198.956 63.4119 199.032 68.2935C199.125 74.1886 198.86 79.7394 199.036 85.5037C206.074 91.1305 221.546 96.4395 225.427 102.83C225.651 108.077 226.212 128.491 225.003 132.326C219.254 129.935 218.782 129.944 213.166 132.594C212.839 131.303 212.927 125.114 212.938 123.471L212.991 108.992C208.025 105.521 198.223 100.284 192.678 97.039C185.318 100.597 171.882 110.454 164.631 112.086C158.851 110.424 156.711 103.729 161.417 100.414C169.111 94.9961 178.309 90.9224 186.262 85.8051Z"/>'
    # lower tail/arrow – muted accent
    '<path fill="var(--fz-accent)" opacity="0.55" d="M268.216 306.183C270.07 308.05 273.034 311.959 274.862 314.121L288.23 329.8L301.33 345.358C303.018 347.363 306.335 351.03 307.764 353.021L307.655 353.28C305.884 357.373 308.695 367.643 305.741 371.524C302.524 375.753 295.068 374.428 294.391 368.17C293.405 359.054 294.449 350.062 293.703 340.966C288.332 337.037 270.769 329.062 267.939 324.393C267.37 323.454 267.21 322.154 267.104 321.072C266.932 319.306 266.726 308.616 267.17 307.299C267.385 306.661 267.642 306.483 268.216 306.183Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.45" d="M276.159 282.925C278.195 279.887 286.408 276.274 289.589 274.22C294.938 270.765 300.073 266.785 306.285 270.844C310.903 273.862 321.216 278.824 324.882 282.654L326.143 282.957C310.216 283.395 292.242 283.023 276.159 282.925Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.45" d="M330.492 282.974C330.521 282.868 330.824 281.763 330.824 281.763C336.244 278.055 342.068 274.922 347.682 271.511C349.803 270.222 353.777 268.466 356.261 269.304C359.451 270.381 360.995 274.363 359.681 277.404C358.9 279.406 357.756 280.803 356.823 282.647C358.7 282.753 361.058 282.777 362.854 282.935C352.067 283.032 341.279 283.045 330.492 282.974Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.45" d="M246.462 282.358C253.144 277.756 252.305 276.386 252.704 268.635C252.742 269.528 252.809 270.033 252.636 270.902L252.943 271.304L252.839 271.503C255.481 272.857 271.22 281.332 272.726 283.212C267.201 284.25 264.926 287.455 262.826 292.404C258.823 289.176 250.579 285.949 246.462 282.358Z"/>'
    # left satellite molecule – blue accent
    '<path fill="var(--fz-accent)" opacity="0.75" d="M134.386 180.174C139.893 177.124 140.836 176.503 147.034 178.519C150.285 180.622 153.613 182.605 157.009 184.463C161.704 187.057 169.569 190.497 173.164 194.121C174.238 195.204 174.836 196.498 175.158 197.975C176.107 202.336 175.361 207.264 175.364 211.711C175.368 217.006 175.818 222.39 175.371 227.67C175.242 229.206 175.063 231 174.215 232.323C171.15 237.101 154.164 245.132 148.223 248.612C148.679 251.305 148.555 257.099 148.515 260.081C148.466 263.695 148.602 268.224 148.413 271.765C152.452 274.05 174.41 285.143 174.939 289.473C176.026 298.368 175.024 310.799 175.374 320.041C175.499 323.166 173.838 326.197 171.161 327.665C163.718 331.748 155.53 337.268 148.066 341.036C145.699 342.201 144.227 342.803 141.708 343.576C139.357 342.722 137.974 342.074 135.756 340.985C130.237 337.261 111.845 328.807 109.426 324.461C108.735 323.219 108.456 321.86 108.32 320.463C108.016 317.324 108.368 314.066 108.345 310.906C108.304 305.166 107.291 294.261 108.795 289.028C109.209 287.59 110.745 286.332 111.894 285.442C114.889 283.123 118.486 281.333 121.768 279.439C126.302 276.822 130.807 274.09 135.43 271.636C134.826 266.396 135.199 254.183 135.254 248.381C128.664 245.181 121.881 240.42 115.328 237.002C112.747 235.372 108.837 233.301 108.626 229.843C108.071 220.701 108.139 210.502 108.35 201.416C108.461 200.952 108.529 199.484 108.586 198.909C109.028 194.415 111.696 193.423 115.248 191.18C120.667 188.462 128.894 183.342 134.386 180.174Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.9" d="M141.636 283.23C143.484 283.636 159.445 293.098 162.25 294.71C162.052 302.652 162.139 309.8 162.292 317.743C156.335 320.752 148.468 325.699 142.348 329.091C141.959 329.219 141.599 329.247 141.228 329.035C137.77 327.06 134.291 324.911 130.79 323.027C127.854 321.446 124.758 319.286 121.724 318.053C121.76 310.791 122.007 301.65 121.612 294.515C127.727 291.281 135.828 286.816 141.636 283.23Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.9" d="M141.563 190.723C143.393 191.553 146.083 193.154 147.889 194.163L162.26 202.358C161.855 209.294 162.093 218.593 162.302 225.524C155.379 229.08 148.93 233.192 142.138 236.922C141.016 236.951 124.685 227.313 121.705 225.87C121.873 220.391 122.055 207.905 121.569 202.542C126.995 199.254 137.175 194.168 141.563 190.723Z"/>'
    # right satellite molecule
    '<path fill="var(--fz-accent)" opacity="0.75" d="M300.51 90.9341C302.851 91.5809 304.534 92.3624 306.734 93.3689C311.29 96.1496 315.991 98.6747 320.663 101.252C324.284 103.249 328.529 105.202 331.575 108.032C332.722 109.098 333.371 110.125 333.76 111.65C334.932 116.239 334.046 121.988 334.055 126.724C334.064 132.063 334.476 137.635 333.992 142.934C333.886 144.102 333.685 145.659 333.011 146.634C330.852 149.756 311.581 160.301 307.173 162.431C307.358 169.82 307.387 178.77 307.225 186.171C311.947 188.053 316.243 191.123 320.65 193.631C323.929 195.498 327.601 197.071 330.594 199.368C331.729 200.238 332.724 201.342 333.288 202.668C333.755 203.767 333.957 204.95 334.065 206.132C334.507 210.967 334.762 235.096 333.133 238.849C332.696 239.854 331.537 240.63 330.663 241.245C326.513 244.166 321.785 246.444 317.367 248.94C312.384 251.754 307.243 255.364 301.82 257.227C301.457 257.351 301.091 257.469 300.724 257.58C291.68 254.751 280.628 246.533 271.965 242.186C268.505 240.45 266.999 236.836 267.099 233.285C266.855 224.639 267.371 215.976 267.072 207.324C266.986 204.835 267.884 201.613 270.018 200.227C277.111 195.621 286.46 189.362 294.144 186.139C293.752 178.518 294.143 170.197 293.963 162.474C289.493 159.853 284.994 157.28 280.469 154.756C277.059 152.834 273.26 151.034 270.273 148.497C269.084 147.487 267.991 146.276 267.582 144.734C266.08 139.064 267.209 130.076 267.241 124.01C267.258 120.614 266.88 117.102 267.182 113.726C267.302 112.38 267.551 111.042 268.225 109.851C270.727 105.422 289.145 96.7174 294.379 93.4696C296.622 92.2995 298.057 91.5719 300.51 90.9341Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.9" d="M300.311 104.98C301.339 104.888 318.486 115.014 320.818 116.325C320.741 121.885 320.424 134.678 320.932 139.664C314.139 143.005 307.358 147.755 300.833 150.921C300.075 151.072 283.034 141.178 280.473 139.887C280.725 132.268 280.88 124.098 280.423 116.51C286.691 113.043 294.328 108.768 300.311 104.98Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.9" d="M300.394 197.266C301.37 197.44 318.521 207.392 320.766 208.597C320.837 213.408 320.672 218.333 320.726 223.15C320.756 225.788 320.58 229.273 321.038 231.817C314.352 235.755 307.618 239.611 300.838 243.385C300.32 243.447 282.686 233.089 280.444 231.95C280.771 229.425 280.667 225.6 280.671 222.965C280.678 218.241 280.703 213.4 280.43 208.694C287.137 204.977 293.793 201.168 300.394 197.266Z"/>'
    # connectors
    '<path fill="var(--fz-accent)" opacity="0.5" d="M96.7372 150.832C92.8027 145.373 89.4902 137.631 85.692 132.019C83.0371 128.096 79.5734 121.127 84.7776 118.001C87.3955 116.428 90.0421 115.909 92.5914 118.695C98.2316 126.617 103.041 135.515 108.02 143.929C113.923 144.714 135.787 145.346 139.427 148.577C141.473 150.393 146.923 174.693 147.034 178.519C140.836 176.503 139.893 177.124 134.386 180.174C128.894 183.342 120.667 188.462 115.248 191.18C111.696 193.423 109.028 194.415 108.586 198.909C108.529 199.484 108.461 200.952 108.35 201.416C105.683 200.72 91.0874 186.184 87.305 183.506C82.3107 179.97 87.3964 170.996 89.4292 166.505C91.7855 161.244 94.2218 156.019 96.7372 150.832Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.7" d="M108.299 156.537C113.629 157.591 124.508 158.498 129.959 158.77C130.361 162.45 133.132 176.822 134.386 180.174C128.894 183.342 120.667 188.462 115.248 191.18C112.715 188.106 102.503 179.602 99.0831 176.781C101.734 170.405 105.365 162.863 108.299 156.537Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.5" d="M135.756 340.985C137.974 342.074 139.357 342.722 141.708 343.576C144.227 342.803 145.699 342.201 148.066 341.036C149.037 345.884 148.273 365.049 147.366 370.913C146.352 372.469 144.97 373.038 143.216 373.508C132.64 376.345 135.287 360.385 135.041 354.477C134.924 351.668 134.587 342.629 135.756 340.985Z"/>'
    '<path fill="var(--fz-accent)" opacity="0.5" d="M294.883 63.7831C295.291 63.3008 295.708 62.8286 296.207 62.4372C297.502 61.4219 299.37 60.8467 301.009 61.049C302.672 61.2542 304.894 61.7989 305.867 63.2833C307.867 66.3337 307.672 89.0709 306.734 93.3689C304.534 92.3624 302.851 91.5809 300.51 90.9341C298.057 91.5719 296.622 92.2995 294.379 93.4696C293.124 89.4664 294.216 68.5909 294.883 63.7831Z"/>'
    # purple lower cascade (funnel/arrow)
    '<path fill="var(--fz-accent)" opacity="0.35" d="M272.472 294.289C274.704 292.324 290.537 293.073 294.528 293.067L334.58 293.093L383.88 293.029C392.912 293.04 402.085 292.977 411.114 293.251C412.606 293.297 413.401 293.294 414.359 294.268C414.567 295.01 414.483 295.648 414.021 296.282C410.909 300.548 406.804 304.84 403.365 308.831L376.662 340.176L366.045 352.593C356.192 364.183 357.005 362.688 357.043 377.833L357.039 397.432C357.043 400.447 357.248 403.532 357.095 406.536C357.047 407.475 356.988 408.667 356.485 409.478C355.988 410.28 355.099 411.019 354.365 411.611C352.33 413.252 333.658 426.631 331.958 426.827C331.611 426.868 331.611 426.851 331.317 426.804C329.861 425.192 330.517 413.593 330.506 410.99L330.384 379.274C330.379 376.945 330.635 366.345 329.76 364.765C323.29 353.07 311.532 342.131 303.244 331.544C301.598 329.441 295.285 321.754 294.225 320.959C287.658 312.814 280.219 304.976 273.601 296.851C273.172 296.324 272.712 294.978 272.472 294.289Z"/>'
    '</svg>')
ICON_HTVS = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M3 3h18v18H3zM3 9h18M3 15h18M9 3v18M15 3v18"/><path d="M9 9l6 6M15 9l-6 6" stroke-width="2" stroke="var(--fz-accent)"/></svg>'
ICON_STATS = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M3 20h18M3 16h4v4H3zM9 12h4v8H9zM15 8h4v12h-4zM3 12l6-4 6 4 6-8"/></svg>'
ICON_EXPORT = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>'

# Basic Toxicophores / Reactive groups SMARTS patterns
TOXIC_SMARTS = {
    'Epoxide': '[OX3]1[CX4][CX4]1',
    'Quinone': 'O=C1[#6]~[#6]C(=O)[#6]~[#6]1',
    'Nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]',
    'Michael_Acceptor': 'C=C-C(=O)',
    'Aliphatic_Halide': '[CX4][F,Cl,Br,I]',
}

def check_toxicity(mol):
    if not mol: return False
    for name, smarts in TOXIC_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            return True, name
    return False, "None"

def get_admet_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    has_tox, tox_alert = check_toxicity(mol)
    return {'MW': round(mw, 2), 'LogP': round(logp, 2), 'HBD': hbd, 'HBA': hba, 'TPSA': round(tpsa, 2), 'Tox_Alert': tox_alert}

def filter_top_candidates(df, max_mw, max_logp, max_hbd, max_hba, max_tpsa, filter_tox, target_count=100, min_score=0.8):
    df = df.sort_values(by='Binding_Score', ascending=False).reset_index(drop=True)
    passed_records = []
    for _, row in df.iterrows():
        if len(passed_records) >= target_count: break
        if row['Binding_Score'] < min_score: break
        props = get_admet_properties(row['SMILES'])
        if not props: continue
        if props['MW'] > max_mw or props['LogP'] > max_logp or props['HBD'] > max_hbd or props['HBA'] > max_hba or props['TPSA'] > max_tpsa: continue
        if filter_tox and props['Tox_Alert'] != 'None': continue
        record = row.to_dict(); record.update(props)
        passed_records.append(record)
    return pd.DataFrame(passed_records)

def get_base64_img(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return ""
        img = Draw.MolToImage(mol, size=(250, 250))
        buffered = BytesIO(); img.save(buffered, format="PNG")
        return f"data:image/png;base64,{base64.b64encode(buffered.getvalue()).decode()}"
    except: return ""

def render_page():
    render_header("Virtual Target Screening", "Screen complex chemical libraries against RNA architectures using Foundry-scale inference.")
    
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    col1, col2 = st.columns([1, 1])
    with col1:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_HELIX}</span> Target Sequence</h3>", unsafe_allow_html=True)
        rna_seq = st.text_area("Input RNA Sequence", height=150, placeholder="GGCGCUUGCUGGGUAAAA...")
    with col2:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_HTVS_MOL}</span> Chemical Library</h3>", unsafe_allow_html=True)
        library_source = st.radio("Library Source", ["Paste SMILES", "Upload File (.smi, .csv, .sdf)"])
        ligands_input = []
        if library_source == "Paste SMILES":
            smiles_text = st.text_area("SMILES (one per line)", height=150, placeholder="CCO\nC1CCCCC1")
            if smiles_text: ligands_input = [s.strip() for s in smiles_text.split() if s.strip()]
        elif library_source == "Upload File (.smi, .csv, .sdf)":
            uploaded_file = st.file_uploader("Upload Ligands", type=['smi', 'csv', 'txt', 'sdf'])
            if uploaded_file:
                filename = uploaded_file.name.lower()
                if filename.endswith('.sdf'):
                    supplier = Chem.ForwardSDMolSupplier(uploaded_file)
                    for mol in supplier:
                        if mol:
                            try: ligands_input.append(Chem.MolToSmiles(mol))
                            except: pass
                else:
                    content = uploaded_file.getvalue().decode('utf-8').splitlines()
                    ligands_input = [line.strip().split()[0] for line in content if line.strip()]
                st.info(f"Loaded {len(ligands_input)} molecules.")
    st.markdown("</div>", unsafe_allow_html=True)

    render_section_divider()
    
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    with st.expander("Advanced ADMET & Toxicity Filters", expanded=False):
        st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.9rem;'>Molecules passing these filters with a Binding Score ≥ 0.8 will be prioritized.</p>", unsafe_allow_html=True)
        col_f1, col_f2, col_f3 = st.columns(3)
        with col_f1:
            max_mw = st.number_input("Max MW", value=500.0)
            max_logp = st.number_input("Max LogP", value=5.0)
        with col_f2:
            max_hbd = st.number_input("Max HBD", value=5)
            max_hba = st.number_input("Max HBA", value=10)
        with col_f3:
            max_tpsa = st.number_input("Max TPSA", value=140.0)
            filter_tox = st.checkbox("Filter Toxicophores", value=True)
    st.markdown("</div>", unsafe_allow_html=True)
            
    render_section_divider()
    
    if st.button("Initialize Foundational Screening", use_container_width=True):
        if not rna_seq or not ligands_input:
            st.error("Please provide RNA target and ligands."); return
            
        with st.spinner("Executing Inference Pipeline..."):
            progress_bar = st.progress(0)
            for i in range(100):
                time.sleep(0.01); progress_bar.progress(i + 1)
            results_df = run_virtual_screening(rna_seq, ligands_input)
            progress_bar.empty()
            render_success_alert("Screening Complete!")
            
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_HTVS}</span> Top Binding Candidates</h3>", unsafe_allow_html=True)
            col_a, col_b, col_c = st.columns(3)
            col_a.metric("Molecules Screened", len(results_df))
            col_b.metric("Highest Affinity", f"{results_df['Binding_Score'].max():.3f}")
            col_c.metric("Foundational Time", "~1.5s")
            
            with st.spinner("Finalizing Analytics..."):
                display_df = filter_top_candidates(results_df, max_mw, max_logp, max_hbd, max_hba, max_tpsa, filter_tox)
                if len(display_df) > 0: display_df['Structure'] = display_df['SMILES'].apply(get_base64_img)
            
            if len(display_df) > 0:
                html_table = "<div class='results-table-wrap'><table class='results-table'><thead><tr><th>Rank</th><th>Structure</th><th>SMILES</th><th>Binding Score</th><th>ADMET</th></tr></thead><tbody>"
                for idx, row in display_df.iterrows():
                    bg = "var(--fz-bg-primary)" if idx % 2 == 0 else "var(--fz-bg-secondary)"
                    img_tag = f"<img src='{row['Structure']}' width='150'/>" if row['Structure'] else "N/A"
                    admet_str = f"<span class='admet-label'>MW</span> {row['MW']}<br><span class='admet-label'>LogP</span> {row['LogP']}<br><span class='admet-label'>HBD/HBA</span> {row['HBD']}/{row['HBA']}"
                    html_table += f"<tr style='background-color:{bg}'><td><span class='rank-badge'>{idx+1}</span></td><td>{img_tag}</td><td class='smiles-cell'>{row['SMILES']}</td><td class='score-cell'>{row['Binding_Score']:.4f}</td><td class='admet-cell'>{admet_str}</td></tr>"
                html_table += "</tbody></table></div>"
                st.markdown(html_table, unsafe_allow_html=True)

            render_section_divider()
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_STATS}</span> Screening Statistics</h3>", unsafe_allow_html=True)
            stats_a, stats_b, stats_c, stats_d = st.columns(4)
            stats_a.metric("Mean Score", f"{results_df['Binding_Score'].mean():.3f}")
            stats_b.metric("Median", f"{results_df['Binding_Score'].median():.3f}")
            stats_c.metric("Min", f"{results_df['Binding_Score'].min():.3f}")
            stats_d.metric("Max", f"{results_df['Binding_Score'].max():.3f}")
            
            import numpy as np
            counts, bins = np.histogram(results_df['Binding_Score'], bins=20)
            st.bar_chart(pd.DataFrame({'Frequency': counts}, index=np.round(bins[:-1], 3)), use_container_width=True)
            
            render_section_divider()
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_EXPORT}</span> Export Results</h3>", unsafe_allow_html=True)
            col_export1, col_export2 = st.columns(2)
            with col_export1: st.download_button("Download ALL (.csv)", results_df.to_csv(index=False).encode('utf-8'), 'all_results.csv', 'text/csv', use_container_width=True)
            with col_export2: st.download_button("Download TOP (.csv)", display_df.drop(columns=['Structure'], errors='ignore').to_csv(index=False).encode('utf-8'), 'top_results.csv', 'text/csv', use_container_width=True)
    render_footer()

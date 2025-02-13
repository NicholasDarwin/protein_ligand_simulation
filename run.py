import http.server
import socketserver
import os

PORT = 8080  # You can change this port if needed

# Path to the HTML file
HTML_FILE_PATH = "index.html"

class MyHTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        # Check if the request path is "/"
        if self.path == "/":
            try:
                # Read the HTML file
                with open(HTML_FILE_PATH, 'r') as f:
                    html_content = f.read()

                # Send response with the HTML content
                self.send_response(200)
                self.send_header('Content-type', 'text/html')
                self.end_headers()
                self.wfile.write(html_content.encode())  # Send the HTML content

            except FileNotFoundError:
                self.send_error(404, "File Not Found")
        else:
            self.send_error(404, "Endpoint not found")

# Start the server
with socketserver.TCPServer(("", PORT), MyHTTPRequestHandler) as httpd:
    print(f"Serving at http://localhost:{PORT}")
    httpd.serve_forever()
